// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointStd_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)

//#include <iostream>

#include <ostream>

#include "../../include/std/EvalNewPointStd_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

ReducedSpaceSQPPack::EvalNewPointStd_Step::EvalNewPointStd_Step(
		const deriv_tester_ptr_t& deriv_tester )
	:
		  deriv_tester_(deriv_tester)
		, new_point_(true)
{}

bool ReducedSpaceSQPPack::EvalNewPointStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;
	using GeneralIterationPack::print_algorithm_step;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPReduced	&nlp	= algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Set the references to the current point's quantities to be updated
	bool f_k_updated = s.f().updated_k(0);
	nlp.set_f(	&s.f().set_k(0)			);
	nlp.set_Gf(	&s.Gf().set_k(0).v()	);
	bool c_k_updated = s.c().updated_k(0);
	nlp.set_c(	&s.c().set_k(0).v()		);
	nlp.set_Gc(	&s.Gc().set_k(0)		);

	VectorWithNorms& x = s.x().get_k(0);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "\n||x||inf = " << x.norm_inf() << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx = \n" << x.v();
	}
	
	// allow multiple updates as defined in NLP, and NLPFirstOrderInfo<...>.
	nlp.set_mult_calc(true);

	// Calculate Gc at x_k
//	std::cout
//		<< "EvaluateNewPoint: x(1) = " << x(1)
//		<< ", new_point_ = " << new_point_ << "\n";
	nlp.calc_Gc( x.v(), new_point_ );

	// Form the decomposition of Gc and update Z and Y
	s.decomp_sys().update_decomp( &s.Gc().get_k(0), &s.Z().set_k(0)
		, &s.Y().set_k(0), &s.U().set_k(0), &s.V().set_k(0) );

	// Calculate the rest of the quantities.  If decomp_sys is a variable
	// reduction decomposition system object, then nlp will be hip to the
	// basis selection and will permute these quantities to that basis.
	// Note that x will already be permuted to the current basis.
	nlp.calc_Gf(x.v(), false);
	if( !c_k_updated )
		nlp.calc_c(x.v(), false);
	if( !f_k_updated )
		nlp.calc_f(x.v(), false);

	// Check the derivatives if we are checking the results
	if( s.check_results() ) {
		
		if( olevel >= PRINT_ALGORITHM_STEPS ) {
			out	<< "\n*** Checking derivatives by finite differences\n";
		}

		const bool result = deriv_tester().finite_diff_check(
			  &nlp, s.Gc().get_k(0), s.Gf().get_k(0)(), x()
			, ( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : 0 );
		if( !result ) {
			throw std::logic_error( "EvalNewPointStd_Step::do_step(...) : "
				"Error, the finite test of the first derivatives of the NLP failed" );
		}
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nf         = "	<< s.f().get_k(0)
			<< "\n||Gf||inf = "	<< s.Gf().get_k(0).norm_inf()
			<< "\n||c||inf  = " << s.c().get_k(0).norm_inf() << "\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		s.Gc().get_k(0).output( out << "\nGc_k = \n" );
		s.Z().get_k(0).output( out << "\nZ_k = \n" );
		s.Y().get_k(0).output( out << "\nY_k = \n" );
		if(algo.nlp().m() > algo.nlp().r()) {
			s.U().get_k(0).output( out << "\nU_k = \n" );
			s.V().get_k(0).output( out << "\nV_k = \n" );
		}
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nGf = \n" << s.Gf().get_k(0)();
		out	<< "\nc = \n" << s.c().get_k(0)() << "\n";
	}

	new_point_ = true;	// Figure this out?
	return true;
}

void ReducedSpaceSQPPack::EvalNewPointStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the new point\n"
		<< L << "Gc_k = Gc(x_k) <: R^n -> R^(n x m)\n"
		<< L << "Reorder Gc and determine con_indep s.t.:\n"
		<< L << "    Gc_k(:,con_indep) <: R^(n x r) has full column rank r\n"
		<< L << "Find:\n"
		<< L << "    Z_k <: R^(n x (n-r))         s.t. Gc_k(:,con_indep)' * Z_k = 0\n"
		<< L << "    Y_k <: R^(n x r)             s.t. [Z_k Y_k] and therefore Gc(:,con_indep)'*Y are nonsigular \n"
		<< L << "    U_k <: R^((m -r) x r)        s.t. U_k = Gc_k(:,con_indep)' * Y_k\n"
		<< L << "    V_k <: R^((m -r) x (n -r))   s.t. V_k = Gc_k(:,con_indep)' * Z_k\n"
		<< L << "Gf_k = Gf(x_k) <: R^n -> R^n\n"
		<< L << "if c_k is not updated c_k = c(x_k) <: R^n -> R^m\n"
		<< L << "if f_k is not updated f_k = f(x_k) <: R^n -> R^1\n"
		<< L << "if check_results == true then check Gc and Gf by finite differences.\n"
		;
}
