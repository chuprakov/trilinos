// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointStd_Step.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_Step.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/include/PermutationOut.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

#ifdef _DEBUG
#include "LinAlgPack/include/PermVecMat.h"
#endif

namespace ReducedSpaceSQPPack {

EvalNewPointStd_Step::EvalNewPointStd_Step(
	const decomp_sys_handler_ptr_t                              &decomp_sys_handler
	,const deriv_tester_ptr_t                                   &deriv_tester
	,const decomp_sys_tester_ptr_t                              &decomp_sys_tester
	,const bounds_tester_ptr_t                                  &bounds_tester
	,EFDDerivTesting                                            fd_deriv_testing
	,DecompositionSystemHandler_Strategy::EDecompSysTesting     decomp_sys_testing
	,DecompositionSystemHandler_Strategy::EDecompSysPrintLevel  decomp_sys_testing_print_level
	)
	:decomp_sys_handler_(decomp_sys_handler)
	,deriv_tester_(deriv_tester)
	,decomp_sys_tester_(decomp_sys_tester)
	,bounds_tester_(bounds_tester)
	,fd_deriv_testing_(fd_deriv_testing)
	,decomp_sys_testing_(decomp_sys_testing)
	,decomp_sys_testing_print_level_(decomp_sys_testing_print_level)
{}

bool EvalNewPointStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using GeneralIterationPack::print_algorithm_step;
	using NLPInterfacePack::NLPFirstOrderInfo;

	rSQPAlgo            &algo   = rsqp_algo(_algo);
	rSQPState           &s      = algo.rsqp_state();
	NLPFirstOrderInfo   &nlp    = dyn_cast<NLPFirstOrderInfo>(algo.nlp());

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if(!nlp.is_initialized())
		nlp.initialize(algo.algo_cntr().check_results());

	const size_type
		n  = nlp.n(),
		nb = nlp.num_bounded_x(),
		m  = nlp.m(),
		mI = nlp.mI();
	size_type
		r  = s.decomp_sys().equ_decomp().size();

	// Get the iteration quantity container objects
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	IterQuantityAccess<value_type>
		&f_iq   = s.f();
	IterQuantityAccess<VectorWithOpMutable>
		&x_iq   = s.x(),
		&nu_iq  = s.nu(),
		*c_iq   = m > 0  ? &s.c() : NULL,
		*h_iq   = mI > 0 ? &s.h() : NULL,
		&Gf_iq  = s.Gf();
	IterQuantityAccess<MatrixWithOp>
		*Gc_iq  = m  > 0 ? &s.Gc() : NULL,
		*Gh_iq  = mI > 0 ? &s.Gh() : NULL,
		*Z_iq   = NULL,
		*Y_iq   = NULL,
		*Uz_iq  = NULL,
		*Uy_iq  = NULL,
		*Vz_iq  = NULL,
		*Vy_iq  = NULL;
	IterQuantityAccess<MatrixWithOpNonsingular>
		*R_iq   = NULL;

	MatrixWithOp::EMatNormType mat_nrm_inf = MatrixWithOp::MAT_NORM_INF;
	
	if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nx is not updated for any k so set x_k = nlp.xinit() ...\n";
		}
		x_iq.set_k(0) = nlp.xinit();
	}
	
	// Validate x
	if( nb && algo.algo_cntr().check_results() ) {
		assert_print_nan_inf(
			x_iq.get_k(0), "x_k", true
			, int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
			);
		if( nlp.num_bounded_x() > 0 ) {
			if(!bounds_tester().check_in_bounds(
				   int(olevel)  >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
				   ,int(olevel) >= int(PRINT_VECTORS)                // print_all_warnings
				   ,int(olevel) >= int(PRINT_ITERATION_QUANTITIES)  // print_vectors
				   ,nlp.xl(),        "xl"
				   ,nlp.xu(),        "xu"
				   ,x_iq.get_k(0),   "x_k"
				   ))
			{
				THROW_EXCEPTION(
					true, TestFailed
					,"EvalNewPointStd_Step::do_step(...) : Error, "
					"the variables bounds xl <= x_k <= xu where violated!" );
			}
		}
	}

	VectorWithOp &x = x_iq.get_k(0);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "\n||x_k||inf = " << x.norm_inf() << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx_k = \n" << x;
	}

	// Set the references to the current point's quantities to be updated
	const bool f_k_updated  = f_iq.updated_k(0);
	const bool Gf_k_updated = Gf_iq.updated_k(0);
	const bool c_k_updated  = m  > 0 ? c_iq->updated_k(0)  : false;
	const bool h_k_updated  = mI > 0 ? h_iq->updated_k(0)  : false;
	const bool Gc_k_updated = m  > 0 ? Gc_iq->updated_k(0) : false;
	const bool Gh_k_updated = mI > 0 ? Gh_iq->updated_k(0) : false;
	nlp.set_f( &f_iq.set_k(0) );
	nlp.set_Gf( &Gf_iq.set_k(0) );
	if( m > 0 ) {
		nlp.set_c( &c_iq->set_k(0) );
		nlp.set_Gc( &Gc_iq->set_k(0) );
	}
	if( mI > 0 ) {
		nlp.set_h( &h_iq->set_k(0) );
		nlp.set_Gh( &Gh_iq->set_k(0) );
	}

	// Allow multiple updates as defined in NLP and NLPFirstOrderInfo interfaces
	nlp.set_multi_calc(true);

	// Calculate Gc and Gh at x_k
	bool new_point = true;
	if(m > 0) {
		nlp.calc_Gc( x, new_point );
		new_point = false;
	}
	if(mI > 0) {
		nlp.calc_Gh( x, new_point );
		new_point = false;
	}

	//
	// Update (or select a new) range/null decomposition
	//
	bool new_decomp_selected = false;
	if( m > 0 ) {
		
		// Update the range/null decomposition
		decomp_sys_handler().update_decomposition(
			algo, s, nlp, decomp_sys_testing(), decomp_sys_testing_print_level()
			,&new_decomp_selected
			);

		r  = s.decomp_sys().equ_decomp().size();

		Z_iq   = ( n > m && r > 0 )      ? &s.Z()  : NULL;
		Y_iq   = ( r > 0 )               ? &s.Y()  : NULL;
		Uz_iq  = ( m  > 0 && m  > r )    ? &s.Uz() : NULL;
		Uy_iq  = ( m  > 0 && m  > r )    ? &s.Uy() : NULL;
		Vz_iq  = ( mI > 0 ) && ( m > 0 ) ? &s.Vz() : NULL;
		Vy_iq  = ( mI > 0 ) && ( m > 0 ) ? &s.Vy() : NULL;
		R_iq   = ( m > 0 )               ? &s.R()  : NULL;

		// Determine if we will test the decomp_sys or not
		const DecompositionSystem::ERunTests
			ds_test_what = ( ( decomp_sys_testing() == DecompositionSystemHandler_Strategy::DST_TEST
							   || ( decomp_sys_testing() == DecompositionSystemHandler_Strategy::DST_DEFAULT
									&& algo.algo_cntr().check_results() ) )
							 ? DecompositionSystem::RUN_TESTS
							 : DecompositionSystem::NO_TESTS );
		
		// Determine the output level for decomp_sys				
		DecompositionSystem::EOutputLevel ds_olevel;
		switch(olevel) {
			case PRINT_NOTHING:
			case PRINT_BASIC_ALGORITHM_INFO:
				ds_olevel = DecompositionSystem::PRINT_NONE;
				break;
			case PRINT_ALGORITHM_STEPS:
			case PRINT_ACTIVE_SET:
				ds_olevel = DecompositionSystem::PRINT_BASIC_INFO;
				break;
			case PRINT_VECTORS:
				ds_olevel = DecompositionSystem::PRINT_VECTORS;
				break;
			case PRINT_ITERATION_QUANTITIES:
				ds_olevel = DecompositionSystem::PRINT_EVERY_THING;
				break;
			default:
				assert(0); // Should not get here!
		};

		// Test the decomposition system
		if( ds_test_what == DecompositionSystem::RUN_TESTS ) {
			// Set the output level
			if( decomp_sys_tester().print_tests() == DecompositionSystemTester::PRINT_NOT_SELECTED ) {
				DecompositionSystemTester::EPrintTestLevel  ds_olevel;
				switch(olevel) {
					case PRINT_NOTHING:
					case PRINT_BASIC_ALGORITHM_INFO:
						ds_olevel = DecompositionSystemTester::PRINT_NONE;
						break;
					case PRINT_ALGORITHM_STEPS:
					case PRINT_ACTIVE_SET:
						ds_olevel = DecompositionSystemTester::PRINT_BASIC;
						break;
					case PRINT_VECTORS:
						ds_olevel = DecompositionSystemTester::PRINT_MORE;
						break;
					case PRINT_ITERATION_QUANTITIES:
						ds_olevel = DecompositionSystemTester::PRINT_ALL;
						break;
					default:
						assert(0); // Should not get here!
				}
				decomp_sys_tester().print_tests(ds_olevel);
				decomp_sys_tester().dump_all( olevel == PRINT_ITERATION_QUANTITIES );
			}
			// Run the tests
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << "\nTesting the range/null decompostion ...\n";
			}
			const bool
				decomp_sys_passed = decomp_sys_tester().test_decomp_system(
					s.decomp_sys()
					,Gc_iq->get_k(0)                   // Gc
					,mI > 0 ? &Gh_iq->get_k(0) : NULL  // Gh
					,&Z_iq->get_k(0)                   // Z
					,&Y_iq->get_k(0)                   // Y
					,&R_iq->get_k(0)                   // R
					,m > r  ? &Uz_iq->get_k(0) : NULL  // Uz
					,m > r  ? &Uy_iq->get_k(0) : NULL  // Uy
					,mI > 0 ? &Vz_iq->get_k(0) : NULL  // Vz
					,mI > 0 ? &Vy_iq->get_k(0) : NULL  // Vy
					,( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : NULL
					);
			THROW_EXCEPTION(
				!decomp_sys_passed, TestFailed
				,"EvalNewPointStd_Step::do_step(...) : Error, "
				"the tests of the decomposition system failed!" );
		}
	}

	// Calculate the rest of the quantities.  If decomp_sys is a variable
	// reduction decomposition system object, then nlp will be hip to the
	// basis selection and will permute these quantities to that basis.
	// Note that x will already be permuted to the current basis.
	nlp.calc_Gf( x, new_point ); new_point = false;
	if( m && (!c_k_updated || new_decomp_selected ) )
		nlp.calc_c( x, false);
	if( mI && (!h_k_updated || new_decomp_selected ) )
		nlp.calc_h( x, false);
	if( !f_k_updated || new_decomp_selected )
		nlp.calc_f( x, false);

	// Check for NaN and Inf
	if(algo.algo_cntr().check_results()) {
		assert_print_nan_inf(f_iq.get_k(0),"f_k",true,&out); 
		if(m)
			assert_print_nan_inf(c_iq->get_k(0),"c_k",true,&out); 
		if(mI)
			assert_print_nan_inf(h_iq->get_k(0),"h_k",true,&out); 
		assert_print_nan_inf(Gf_iq.get_k(0),"Gf_k",true,&out); 
	}

	// Print the iteration quantities before we test the derivatives for debugging

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "\nPrinting the updated iteration quantities ...\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nf_k         = "     << f_iq.get_k(0);
		out << "\n||Gf_k||inf = "     << Gf_iq.get_k(0).norm_inf();
		if(m) {
			out << "\n||c_k||inf  = " << c_iq->get_k(0).norm_inf();
			out << "\n||Gc_k||inf = " << Gc_iq->get_k(0).calc_norm(mat_nrm_inf).value;
			out << "\n||Z||inf    = " << Z_iq->get_k(0).calc_norm(mat_nrm_inf).value;
			out << "\n||Y||inf    = " << Y_iq->get_k(0).calc_norm(mat_nrm_inf).value;
			out << "\n||R||inf    = " << R_iq->get_k(0).calc_norm(mat_nrm_inf).value;
			if(algo.algo_cntr().calc_conditioning()) {
				out << "\ncond_inf(R) = " << R_iq->get_k(0).calc_cond_num(mat_nrm_inf).value;
			}
			if( m > r ) {
				out << "\n||Uz_k||inf = " << Uz_iq->get_k(0).calc_norm(mat_nrm_inf).value;
				out << "\n||Uy_k||inf = " << Uy_iq->get_k(0).calc_norm(mat_nrm_inf).value;
			}
		}
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		out << "\nGc_k =\n" << Gc_iq->get_k(0);
		if(mI)
			out << "\nGh_k =\n" << Gh_iq->get_k(0);
		out << "\nZ_k =\n" << Z_iq->get_k(0);
		out << "\nY_k =\n" << Y_iq->get_k(0);
		out << "\nR_k =\n" << R_iq->get_k(0);
		if( m > r ) {
			out << "\nUz_k =\n" << Uz_iq->get_k(0);
			out << "\nUy_k =\n" << Uy_iq->get_k(0);
		}
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nGf_k = \n" << Gf_iq.get_k(0);
		if(m)
			out	<< "\nc_k = \n" << c_iq->get_k(0);
		out << std::endl;
	}

	// Check the derivatives if we are checking the results
	if(		fd_deriv_testing() == FD_TEST
		|| ( fd_deriv_testing() == FD_DEFAULT && algo.algo_cntr().check_results() )  )
	{
		
		if( olevel >= PRINT_ALGORITHM_STEPS ) {
			out	<< "\n*** Checking derivatives by finite differences\n";
		}

		const bool
			nlp_passed = deriv_tester().finite_diff_check(
				&nlp
				,x
				,nb ? &nlp.xl() : NULL
				,nb ? &nlp.xu() : NULL
				,m  ? &Gc_iq->get_k(0) : NULL
				,mI ? &Gh_iq->get_k(0) : NULL
				,&Gf_iq.get_k(0)
				,olevel >= PRINT_VECTORS
				,( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : NULL
				);
		THROW_EXCEPTION(
			!nlp_passed, TestFailed
			,"EvalNewPointStd_Step::do_step(...) : Error, "
			"the tests of the nlp derivatives failed!" );
	}

	return true;
}

void EvalNewPointStd_Step::print_step(
	 const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	const rSQPAlgo   &algo = rsqp_algo(_algo);
	const rSQPState  &s    = algo.rsqp_state();
	out
		<< L << "*** Evaluate the new point and update the range/null decomposition\n"
		<< L << "if nlp is not initialized then initialize the nlp\n"
		<< L << "if x is not updated for any k then set x_k = xinit\n"
		<< L << "if m > 0 and Gc_k is not updated Gc_k = Gc(x_k) <: space_x|space_c\n"
		<< L << "if mI > 0 Gh_k is not updated Gh_k = Gh(x_k) <: space_x|space_h\n"
		<< L << "if m > 0 then\n"
		<< L << "  For Gc_k = [ Gc_k(:,equ_decomp), Gc_k(:,equ_undecomp) ] where:\n"
		<< L << "    Gc_k(:,equ_decomp) <: space_x|space_c(equ_decomp) has full column rank r\n"
		<< L << "  Find:\n"
		<< L << "    Z_k  <: space_x|space_null    s.t. Gc_k(:,equ_decomp)' * Z_k = 0\n"
		<< L << "    Y_k  <: space_x|space_range   s.t. [Z_k Y_k] is nonsigular \n"
		<< L << "    R_k  <: space_c(equ_decomp)|space_range\n"
		<< L << "                                  s.t. R_k = Gc_k(:,equ_decomp)' * Y_k\n"
		<< L << "    if m > r : Uz_k <: space_c(equ_undecomp)|space_null\n"
		<< L << "                                  s.t. Uz_k = Gc_k(:,equ_undecomp)' * Z_k\n"
		<< L << "    if m > r : Uy_k <: space_c(equ_undecomp)|space_range\n"
		<< L << "                                  s.t. Uy_k = Gc_k(:,equ_undecomp)' * Y_k\n"
		<< L << "    if mI > 0 : Vz_k <: space_h|space_null\n"
		<< L << "                                  s.t. Vz_k = Gh_k' * Z_k\n"
		<< L << "    if mI > 0 : Vy_k <: space_h|space_range\n"
		<< L << "                                  s.t. Vy_k = Gh_k' * Y_k\n"
		<< L << "  begin update decomposition (class \'" << typeid(decomp_sys_handler()).name() << "\')\n"
		;
	decomp_sys_handler().print_update_decomposition( algo, s, out, L + "    " );
	out
		<< L << "  end update decomposition\n"
		<< L << "  if ( (decomp_sys_testing==DST_TEST)\n"
		<< L << "    or (decomp_sys_testing==DST_DEFAULT and check_results==true)\n"
		<< L << "    ) then\n"
		<< L << "    check properties for Z_k, Y_k, R_k, Uz_k, Uy_k, Vz_k and Vy_k.\n"
		<< L << "  end\n"
		<< L << "end\n"
		<< L << "Gf_k = Gf(x_k) <: space_x\n"
		<< L << "if m > 0 and c_k is not updated c_k = c(x_k) <: space_c\n"
		<< L << "if mI > 0 and h_k is not updated h_k = h(x_k) <: space_h\n"
		<< L << "if f_k is not updated f_k = f(x_k) <: REAL\n"
		<< L << "if ( (fd_deriv_testing==FD_TEST)\n"
		<< L << "  or (fd_deriv_testing==FD_DEFAULT and check_results==true)\n"
		<< L << "  ) then\n"
		<< L << "  check Gc_k (if m > 0), Gh_k (if mI > 0) and Gf_k by finite differences.\n"
		<< L << "end\n"
		;
}

} // end namespace ReducedSpaceSQPPack
