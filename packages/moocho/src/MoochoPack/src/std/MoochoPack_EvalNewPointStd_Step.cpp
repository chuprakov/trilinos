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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)

#include <ostream>

#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_Step.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace ReducedSpaceSQPPack {

EvalNewPointStd_Step::EvalNewPointStd_Step(
	const deriv_tester_ptr_t&         deriv_tester
	,const decomp_sys_tester_ptr_t&   decomp_sys_tester
	,const bounds_tester_ptr_t&       bounds_tester
	,EFDDerivTesting                  fd_deriv_testing
	,EDecompSysTesting                decomp_sys_testing
	)
	:deriv_tester_(deriv_tester)
	,decomp_sys_tester_(decomp_sys_tester)
	,bounds_tester_(bounds_tester)
	,fd_deriv_testing_(fd_deriv_testing)
	,decomp_sys_testing_(decomp_sys_testing)
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
		nlp.initialize();

	const size_type
		n  = nlp.n(),
		nb = nlp.num_bounded_x(),
		m  = nlp.m(),
		mI = nlp.mI(),
		r  = s.decomp_sys().con_decomp().size();

	IterQuantityAccess<value_type>
		&f_iq   = s.f();
	IterQuantityAccess<VectorWithOpMutable>
		&x_iq   = s.x(),
		*c_iq   = m > 0  ? &s.c() : NULL,
		*h_iq   = mI > 0 ? &s.h() : NULL,
		&Gf_iq  = s.Gf();
	IterQuantityAccess<MatrixWithOp>
		*Gc_iq  = m  > 0           ? &s.Gc() : NULL,
		*Gh_iq  = mI > 0           ? &s.Gh() : NULL,
		*Z_iq   = m  > 0           ? &s.Z()  : NULL,
		*Y_iq   = m  > 0           ? &s.Y()  : NULL,
		*Uz_iq  = m  > 0 && m  > r ? &s.Uz() : NULL,
		*Uy_iq  = m  > 0 && m  > r ? &s.Uy() : NULL,
		*Vz_iq  = m  > 0 && mI > 0 ? &s.Vz() : NULL,
		*Vy_iq  = m  > 0 && mI > 0 ? &s.Vy() : NULL;
	IterQuantityAccess<MatrixWithOpNonsingular>
		*R_iq  = m   > 0           ? &s.R()  : NULL;

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
		out << "\n||x||inf = " << x.norm_inf() << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx = \n" << x;
	}

	// Set the references to the current point's quantities to be updated
	const bool f_k_updated  = f_iq.updated_k(0);
	const bool Gf_k_updated = Gf_iq.updated_k(0);
	const bool c_k_updated  = m  > 0 ? c_iq->updated_k(0)  : false;
	const bool h_k_updated  = mI > 0 ? h_iq->updated_k(0)  : false;
	const bool Gc_k_updated = m  > 0 ? Gc_iq->updated_k(0) : false;
	const bool Gh_k_updated = mI > 0 ? Gh_iq->updated_k(0) : false;
	if(f_k_updated)
		nlp.set_f( NULL );
	else
		nlp.set_f( &f_iq.set_k(0) );
	if(Gf_k_updated)
		nlp.set_Gf( NULL );
	else
		nlp.set_Gf( &Gf_iq.set_k(0) );
	if( m > 0 ) {
		if( c_k_updated )
			nlp.set_c( NULL );
		else
			nlp.set_c( &c_iq->set_k(0) );
		if( Gc_k_updated )
			nlp.set_Gc( NULL );
		else
			nlp.set_Gc( &Gc_iq->set_k(0) );
	}
	if( mI > 0 ) {
		if( h_k_updated )
			nlp.set_h( NULL );
		else
			nlp.set_h( &h_iq->set_k(0) );
		if( Gh_k_updated )
			nlp.set_Gh( NULL );
		else
			nlp.set_Gh( &Gh_iq->set_k(0) );
	}

	// allow multiple updates as defined in NLP and NLPFirstOrderInfo
	nlp.set_mult_calc(true);

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

	if( m > 0 ) {
		
		// ToDo: Deal with variable and constraint permutations?

		// Determine if we will test the decomp_sys or not
		const DecompositionSystem::ERunTests
			 ds_test_what = ( ( decomp_sys_testing() == DST_TEST
							|| ( decomp_sys_testing() == DST_DEFAULT
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

		// Form the decomposition of Gc and Gh and update the decomposition system matrices
		s.decomp_sys().update_decomp(
			&out                               // out
			,ds_olevel                         // olevel
			,ds_test_what                      // test_what
			,Gc_iq->get_k(0)                   // Gc
			,mI > 0 ? &Gh_iq->get_k(0) : NULL  // Gh
			,&Z_iq->set_k(0)                   // Z
			,&Y_iq->set_k(0)                   // Y
			,&R_iq->set_k(0)                   // R
			,m > r  ? &Uz_iq->set_k(0) : NULL  // Uz
			,m > r  ? &Uy_iq->set_k(0) : NULL  // Uy
			,mI > 0 ? &Vz_iq->set_k(0) : NULL  // Vz
			,mI > 0 ? &Vy_iq->set_k(0) : NULL  // Vy
			,DecompositionSystem::MATRICES_ALLOW_DEP_IMPS // ToDo: Change this to MATRICES_INDEP_IMPS
			);

		// Print the matrices before we test for debugging
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
			out << "\nGc_k =\n" << Gc_iq->get_k(0);
			if(mI)
				out << "\nGh_k =\n" << Gh_iq->get_k(0);
			out << "\nZ_k =\n" << Z_iq->get_k(0);
			out << "\nY_k =\n" << Y_iq->get_k(0);
			if( m > r ) {
				out << "\nUz_k =\n" << Uz_iq->get_k(0);
				out << "\nUy_k =\n" << Uy_iq->get_k(0);
			}
			if( mI > r ) {
				out << "\nVz_k =\n" << Vz_iq->get_k(0);
				out << "\nVy_k =\n" << Vy_iq->get_k(0);
			}
		}

		// Test the decomposition system
		if( ds_test_what == DecompositionSystem::RUN_TESTS ) {

			// ToDo: set the print levels if not already set!

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
	if( m && !c_k_updated )
		nlp.calc_c( x, false);
	if( mI && !h_k_updated )
		nlp.calc_h( x, false);
	if( !f_k_updated )
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
		out	<< "\nf_k         = "     << f_iq.get_k(0);
		out << "\n||Gf_k||inf = "     << Gf_iq.get_k(0).norm_inf();
		if(m)
			out << "\n||c_k||inf  = " << c_iq->get_k(0).norm_inf();
		if(mI)
			out << "\n||h_k||inf  = " << c_iq->get_k(0).norm_inf();
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nGf = \n" << Gf_iq.get_k(0);
		if(m)
			out	<< "\nc = \n" << c_iq->get_k(0);
		if(mI)
			out	<< "\nh = \n" << h_iq->get_k(0);
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
				,algo.algo_cntr().max_var_bounds_viol()
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
		<< L << "  *** ToDo: Work on variable and constraint permutations!\n"
		<< L << "  For Gc_k = [ Gc_k(:,con_decomp), Gc_k(:,con_undecomp) ] where:\n"
		<< L << "    Gc_k(:,con_decomp) <: space_x|space_c(con_decomp) has full column rank r\n"
		<< L << "  Find:\n"
		<< L << "    Z_k  <: space_x|space_null    s.t. Gc_k(:,con_decomp)' * Z_k = 0\n"
		<< L << "    Y_k  <: space_x|space_range   s.t. [Z_k Y_k] is nonsigular \n"
		<< L << "    R_k  <: space_c(con_decomp)|space_range\n"
		<< L << "                                  s.t. R_k = Gc_k(:,con_decomp)' * Y_k\n"
		<< L << "    if m > r : Uz_k <: space_c(con_undecomp)|space_null\n"
		<< L << "                                  s.t. Uz_k = Gc_k(:,con_undecomp)' * Z_k\n"
		<< L << "    if m > r : Uy_k <: space_c(con_undecomp)|space_range\n"
		<< L << "                                  s.t. Uy_k = Gc_k(:,con_undecomp)' * Y_k\n"
		<< L << "    if mI > 0 : Vz_k <: space_h|space_null\n"
		<< L << "                                  s.t. Vz_k = Gh_k' * Z_k\n"
		<< L << "    if mI > 0 : Vy_k <: space_h|space_range\n"
		<< L << "                                  s.t. Vy_k = Gh_k' * Y_k\n"
		<< L << "  begin update decomposition (class \'" << typeid(s.decomp_sys()).name() << "\')\n"
		;
	s.decomp_sys().print_update_decomp( out, L + "    " );
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
