// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStep_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/LineSearchFullStep_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
	#include "LinAlgPack/include/assert_print_nan_inf.h"

bool ReducedSpaceSQPPack::LineSearchFullStep_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;
	using LinAlgPack::V_VpV;
	using LinAlgPack::assert_print_nan_inf;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLP			&nlp	= algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "f_k        = " << s.f().get_k(0) << std::endl
			<< "||c_k||inf = " << s.c().get_k(0).norm_inf() << std::endl;
	}

	// alpha_k = 1.0
	if( !s.alpha().updated_k(0) )
		s.alpha().set_k(0) = 1.0;
	
	// x_kp1 = x_k + d_k
	Vector &x_kp1 = s.x().set_k(+1).v();
	x_kp1 = s.x().get_k(0)();
	LinAlgPack::Vp_StV( &x_kp1(), s.alpha().get_k(0), s.d().get_k(0)() );

	// Calcuate f and c at the new point.
	nlp.set_c( &s.c().set_k(+1).v() );
	nlp.calc_c( x_kp1 );
	nlp.set_f( &s.f().set_k(+1) );
	nlp.calc_f( x_kp1, false );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "alpha_k      = " << s.alpha().get_k(0) << std::endl
			<< "||x_kp1||inf = " << s.x().get_k(+1).norm_inf() << std::endl
			<< "f_kp1        = " << s.f().get_k(+1) << std::endl
			<< "||c_kp1||inf = " << s.c().get_k(+1).norm_inf() << std::endl;
	}

	assert_print_nan_inf( s.f().get_k(+1), "f(x_kp1)", true, &out );
	assert_print_nan_inf( s.c().get_k(+1)(), "c(x_kp1)", true, &out );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx_kp1 =\n" << s.x().get_k(+1)()
			<< "\nc_kp1 =\n" << s.c().get_k(+1)(); 
	}

	return true;
}

void ReducedSpaceSQPPack::LineSearchFullStep_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if alpha_k is not updated then\n"
		<< L << "    alpha_k = 1.0\n"
		<< L << "end\n"
		<< L << "x_kp1 = x_k + alpha_k * d_k\n"
		<< L << "f_kp1 = f(x_kp1)\n"
		<< L << "c_kp1 = c(x_kp1)\n";
}
