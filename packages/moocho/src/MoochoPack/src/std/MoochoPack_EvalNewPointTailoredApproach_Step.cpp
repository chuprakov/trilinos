// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_Step.cpp
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

//#include <iostream>

#include <ostream>

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/MatrixIdentConcatStd.h"
#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::EvalNewPointTailoredApproach_Step(
	const deriv_tester_ptr_t      &deriv_tester
	,const bounds_tester_ptr_t    &bounds_tester
	, EFDDerivTesting             fd_deriv_testing
	)
	:deriv_tester_(deriv_tester)
	,bounds_tester_(bounds_tester)
	,fd_deriv_testing_(fd_deriv_testing)
{}

bool ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_MtV;
	using GeneralIterationPack::print_algorithm_step;
	namespace rcp = ReferenceCountingPack;

	rSQPAlgo             &algo   = rsqp_algo(_algo);
	rSQPState            &s      = algo.rsqp_state();
	NLPFirstOrderDirect  &nlp    = dyn_cast<NLPFirstOrderDirect>(algo.nlp());

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
		m  = nlp.m(),
		mI = nlp.mI(),
		r  = nlp.var_dep().size();

	IterQuantityAccess<VectorWithOpMutable>
		&x_iq = s.x();

	if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nx is not updated for any k so set x_k = nlp.xinit() ...\n";
		}
		x_iq.set_k(0) = nlp.xinit();
	}

	// ToDo: Incorporate undecomposed equality constraints in the futrue!

	// Validate x
	if(algo.algo_cntr().check_results()) {
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
					,"EvalNewPointTailoredApproach_Step::do_step(...) : Error, "
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

	// allow multiple updates as defined in NLP
	nlp.set_mult_calc(true);

	// If c_k is not updated then we must compute it
	bool recalc_c = true;
	
	if( !s.c().updated_k(0) ) {
		s.c().set_k(0);
		recalc_c = true;
	}
	else {
		recalc_c = false;
	}
		
	// Get references to Z and Y
	MatrixWithOp
		&Z_k = s.Z().set_k(0),
		&Y_k = s.Y().set_k(0);
	MatrixIdentConcatStd
		&cZ_k = dyn_cast<MatrixIdentConcatStd>(Z_k),
		&cY_k = dyn_cast<MatrixIdentConcatStd>(Y_k);
	cY_k.set_uninitialized(); // Release reference to D in case Y references it.
	// If Z has not been initialized or Z.D is being shared by someone else we need to reconstruct Z.D
	bool reconstruct_Z_D = (cZ_k.rows() == n || cZ_k.cols() != n-r || cZ_k.D_ptr().count() > 1);
	MatrixIdentConcatStd::D_ptr_t
		D_ptr = NULL;
	if( reconstruct_Z_D )
		D_ptr = nlp.space_D()->create_member();
	else
		D_ptr = cZ_k.D_ptr();

	// Compute all the quantities.
	VectorWithOpMutable
		&py_k  = s.py().set_k(0);
	nlp.calc_point(
		x
		,!s.f().updated_k(0) ? &s.f().set_k(0) : (value_type*)NULL
		,&s.c().get_k(0)
		,recalc_c
		,(mI && !s.h().updated_k(0)) ? &s.h().set_k(0) : (VectorWithOpMutable*)NULL 
		,&s.Gf().set_k(0)
		,&py_k                                   // -inv(C)*c
		,&s.rGf().set_k(0)
		,NULL                                    // GcU
		,mI ? &s.Gh().set_k(0) : NULL
		,const_cast<MatrixWithOp*>(D_ptr.get())  // -inv(C)*N
		,NULL                                    // V
		,NULL                                    // P
		);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "\nf_k           = " << s.f().get_k(0);
		out << "\n||c_k||inf    = " << s.c().get_k(0).norm_inf();
		if(mI)
			out << "\n||h_k||inf    = " << s.h().get_k(0).norm_inf();
		out << "\n||Gf_k||inf   = " << s.Gf().get_k(0).norm_inf();
		out << "\n||py_k||inf   = " << s.py().get_k(0).norm_inf();
		out << "\n||rGf_k||inf  = " << s.rGf().get_k(0).norm_inf();
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nc_k  = \n" << s.c().get_k(0);
		if(mI)
			out	<< "\nh_k  = \n" << s.h().get_k(0);
		out	<< "\nGf_k = \n" << s.Gf().get_k(0);
		out	<< "\npy_k = \n" << s.py().get_k(0);
		out	<< "\nrGf_k = \n" << s.rGf().get_k(0);
		out << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		out << "\nD = -inv(C)*N = \n" << *D_ptr;
		out << std::endl;
	}

	if(algo.algo_cntr().check_results()) {
		assert_print_nan_inf(s.f().get_k(0),   "f_k",true,&out); 
		assert_print_nan_inf(s.c().get_k(0),   "c_k",true,&out); 
		if(mI)
			assert_print_nan_inf(s.c().get_k(0),"h_k",true,&out); 
		assert_print_nan_inf(s.Gf().get_k(0),  "Gf_k",true,&out); 
		assert_print_nan_inf(s.py().get_k(0),  "py_k",true,&out);
		assert_print_nan_inf(s.rGf().get_k(0), "rGf_k",true,&out);
	}

	// Check the derivatives if we are checking the results
	if(		fd_deriv_testing() == FD_TEST
		|| ( fd_deriv_testing() == FD_DEFAULT && algo.algo_cntr().check_results() )  )
	{
		
		if( olevel >= PRINT_ALGORITHM_STEPS ) {
			out	<< "\n*** Checking derivatives by finite differences\n";
		}
		
		const bool has_bounds = nlp.num_bounded_x() > 0;
		const bool result = deriv_tester().finite_diff_check(
			&nlp
			,x
			,has_bounds ? &nlp.xl() : (const VectorWithOp*)NULL
			,has_bounds ? &nlp.xu() : (const VectorWithOp*)NULL
			,algo.algo_cntr().max_var_bounds_viol()
			,&s.c().get_k(0)
			,mI ? &s.h().get_k(0) : (const VectorWithOp*)NULL
			,&s.Gf().get_k(0)
			,&s.py().get_k(0)
			,&s.rGf().get_k(0)
			,NULL                                                // GcU
			,mI ? &s.Gh().get_k(0) : (const MatrixWithOp*)NULL
			,D_ptr.get()
			,NULL                                                // V
			,NULL                                                // P
			,olevel >= PRINT_VECTORS
			,( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : (std::ostream*)NULL
			);
		if( !result ) {
			throw std::logic_error( "EvalNewPointTailoredApproach_Step::do_step(...) : "
				"Error, the finite derivative test of the first derivatives of the NLP failed" );
		}
	}

	if( reconstruct_Z_D ) {
		//
		// Z = [      D     ] space_xD
		//     [      I     ] space_xI
		//        space_xI
		//
		cZ_k.initialize(
			nlp.space_x()                                          // space_cols
			,nlp.space_x()->sub_space(nlp.var_indep())->clone()    // space_rows
			,MatrixIdentConcatStd::TOP                             // top_or_bottom
			,1.0                                                   // alpha
			,D_ptr                                                 // D_ptr
			,BLAS_Cpp::no_trans                                    // D_trans
			);
	}

	// Compute py and Y
	calc_py_Y( nlp, D_ptr, &py_k, &cY_k, olevel, out ); 

	// Compute Ypy = Y*py
	VectorWithOpMutable
		&Ypy_k  = s.Ypy().set_k(0);
	V_MtV( &Ypy_k, cY_k, BLAS_Cpp::no_trans, py_k );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||py_k||inf   = " << s.py().get_k(0).norm_inf();
		out << "\n||Ypy_k||inf  = " << s.Ypy().get_k(0).norm_inf();
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\npy_k = \n"  << s.py().get_k(0);
		out	<< "\nYpy_k = \n" << s.Ypy().get_k(0);
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		out << "\nZ_k = \n" << s.Z().get_k(0);
		out << "\nY_k = \n" << s.Y().get_k(0);
		out << std::endl;
	}

	return true;

}

void ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the new point for the \"Tailored Approach\"\n"
		<< L << "if nlp is not initialized then initialize the nlp\n"
		<< L << "if x is not updated for any k then set x_k = xinit\n"
		<< L << "if f_k is not updated f_k = f(x_k) <: R^n -> R^1\n"
		<< L << "if c_k is not updated c_k = c(x_k) <: R^n -> R^m\n"
		<< L << "if mI > 0 and h_k is not updated h_k = h(x_k) <: R^n -> R^mI\n"
		<< L << "Gf_k = Gf(x_k) <: R^n -> R^n\n"
		<< L << "For Gc = [ C' ; N' ] = Gc(x_k) <: R^n -> R^(n x m) compute:\n"
		<< L << "    py_k = -inv(C)*c_k\n"
		<< L << "    D = -inv(C)*N <: R^(n x (n-m))\n"
		<< L << "    rGf_k = Gf_k(var_indep) + D'*Gf_k(var_dep)\n"
		<< L << "    Z_k = [ D ; I ] <: R^(n x (n-m))\n"
		<< L << "if (fd_deriv_testing==FD_TEST) or (fd_deriv_testing==FD_DEFAULT and check_results==true) then\n"
		<< L << "    check Gf_k, py_k, rGf_k, and D by finite differences.\n"
		<< L << "end\n";
	print_calc_py_Y( out, L );
	out
		<< L << "Ypy_k = Y_k * py_k\n";
}
