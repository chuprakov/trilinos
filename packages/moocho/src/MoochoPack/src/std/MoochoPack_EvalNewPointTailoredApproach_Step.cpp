// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)

//#include <iostream>

#include <ostream>

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"
#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproach.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "ConstrainedOptimizationPack/include/DenseIdentVertConcatMatrixSubclass.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::EvalNewPointTailoredApproach_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, EFDDerivTesting				fd_deriv_testing
		)
	:
		  deriv_tester_(deriv_tester)
		, fd_deriv_testing_(fd_deriv_testing)
{}

bool ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgPack::norm_inf;
	using LinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_MtV;
	using GeneralIterationPack::print_algorithm_step;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPrSQPTailoredApproach
				&nlp	= dyn_cast<NLPrSQPTailoredApproach>(algo.nlp());

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// ToDo: Incorroparate undecomposed equality constraints in the futrue!

	// ToDo: Put this under check results?
	assert_print_nan_inf(s.x().get_k(0)(), "x_k",true,&out); 

	VectorWithNorms& x = s.x().get_k(0);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "\n||x||inf = " << x.norm_inf() << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx = \n" << x.v();
	}

	// allow multiple updates as defined in NLP, and NLPFirstOrderInfo<...>.
	nlp.set_mult_calc(true);

	// If c_k is not updated then we must set size to zero so it will be
	// computed
	if( !s.c().updated_k(0) )
		s.c().set_k(0).v().resize(0);

	// Get a reference to D = -inv(C)*N storage in Z = [ D; I ].
	MatrixWithOp
		&Z_k = s.Z().set_k(0);
	DenseIdentVertConcatMatrixSubclass
		&cZ_k = dyn_cast<DenseIdentVertConcatMatrixSubclass>(Z_k);
	cZ_k.m().resize( nlp.n(), nlp.n() - nlp.r(), true);	// Z = [ D; I ]
	GenMatrixSlice
		D = cZ_k.m().D();

	// Compute all the quantities.
	nlp.calc_point(
		  x()
		, s.f().updated_k(0) ? (value_type*)NULL : &s.f().set_k(0)
		, &s.c().get_k(0).v() 	// If we need to compute c then c.size() == 0
		, &s.Gf().set_k(0).v()
		, &s.py().set_k(0).v()	// -inv(C)*c
		, &D					// -inv(C)*N
		);

	// ToDo: Put this under check results?
	assert_print_nan_inf(s.f().get_k(0), "f_k",true,&out); 
	assert_print_nan_inf(s.c().get_k(0)(), "c_k",true,&out); 
	assert_print_nan_inf(s.Gf().get_k(0)(), "Gf_k",true,&out); 
	assert_print_nan_inf(s.py().get_k(0)(), "py_k",true,&out); 

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out
			<< "\nf_k           = "	<< s.f().get_k(0)
			<< "\n||Gf_k||inf   = "	<< s.Gf().get_k(0).norm_inf()
			<< "\n||c_k||inf    = " << s.c().get_k(0).norm_inf()
			<< "\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nGf_k = \n" << s.Gf().get_k(0)();
		out	<< "\nc_k  = \n" << s.c().get_k(0)() << "\n";
	}

	// Check the derivatives if we are checking the results
	if(		fd_deriv_testing() == FD_TEST
		|| ( fd_deriv_testing() == FD_DEFAULT && algo.algo_cntr().check_results() )  )
	{
		
		if( olevel >= PRINT_ALGORITHM_STEPS ) {
			out	<< "\n*** Checking derivatives by finite differences\n";
		}

		const bool result = deriv_tester().finite_diff_check(
			  &nlp, x(), &nlp.xl(), &nlp.xu(), algo.algo_cntr().max_var_bounds_viol()
			, &D, &s.py().get_k(0)(), &s.Gf().get_k(0)(), &s.c().get_k(0)()
			, olevel >= PRINT_VECTORS
			, ( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : 0 );
		if( !result ) {
			throw std::logic_error( "EvalNewPointTailoredApproach_Step::do_step(...) : "
				"Error, the finite derivative test of the first derivatives of the NLP failed" );
		}
	}

	// Compute py and Ypy
	calc_py_Ypy( D, &s.py().get_k(0).v()(), &s.Ypy().set_k(0).v(), olevel, out ); 

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||py_k||inf   = " << s.py().get_k(0).norm_inf()
			<< "\n||Ypy_k||inf  = " << s.py().get_k(0).norm_inf()
			<< "\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\npy_k = \n" << s.py().get_k(0)();
		out	<< "\nYpy_k = \n" << s.Ypy().get_k(0)()	<< "\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		s.Z().get_k(0).output( out << "\nZ_k = \n" );
	}

	// For now we just need to update this to nothing.
	s.V().set_k(0);

	return true;

}

void ReducedSpaceSQPPack::EvalNewPointTailoredApproach_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the new point for the \"Tailored Approach\"\n"
		<< L << "if c_k is not updated c_k = c(x_k) <: R^n -> R^m\n"
		<< L << "if f_k is not updated f_k = f(x_k) <: R^n -> R^1\n"
		<< L << "Gf_k = Gf(x_k) <: R^n -> R^n\n"
		<< L << "For Gc = [ C' ; N' ] = Gc(x_k) <: R^n -> R^(n x m) compute:\n"
		<< L << "    py_k = -inv(C)*c_k\n"
		<< L << "    D = -inv(C)*N <: R^(n x (n-m))\n"
		<< L << "    Z_k = [ D ; I ] <: R^(n x (n-m))\n"
		<< L << "if (fd_deriv_testing==FD_TEST) or (fd_deriv_testing==FD_DEFAULT and check_results==true) then\n"
		<< L << "    check Gf_k, py_k, and D by finite differences.\n"
		<< L << "end\n";
	print_calc_Y_py_Ypy( out, L );
}
