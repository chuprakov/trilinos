// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStep_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/LineSearchFullStep_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

bool ReducedSpaceSQPPack::LineSearchFullStep_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;
	using LinAlgPack::V_VpV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPReduced	&nlp	= algo.nlp();

	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "f_k        = " << s.f().get_k(0) << std::endl
			<< "||c_k||inf = " << norm_inf( s.c().get_k(0) ) << std::endl;
	}

	// alpha_k = 1.0
	s.alpha().set_k(0) = 1.0;
	
	// x_kp1 = x_k + d_k
	Vector &x_kp1 = s.x().set_k(+1);
	LinAlgPack::V_VpV( &x_kp1, s.x().get_k(0), s.d().get_k(0)() );

	// Calcuate f and c at the new point.
	nlp.set_c( &s.c().set_k(+1) );
	nlp.calc_c( x_kp1 );
	nlp.set_f( &s.f().set_k(+1) );
	nlp.calc_f( x_kp1, false );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "alpha_k      = " << s.f().get_k(0) << std::endl
			<< "||x_kp1||inf = " << norm_inf( s.x().get_k(+1) ) << std::endl
			<< "f_kp1        = " << s.f().get_k(+1) << std::endl
			<< "||c_kp1||inf = " << norm_inf( s.c().get_k(+1) ) << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx_kp1 =\n" << s.x().get_k(+1)
			<< "\nc_kp1 =\n" << s.c().get_k(+1); 
	}

	return true;
}

void ReducedSpaceSQPPack::LineSearchFullStep_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "alpha_k = 1.0\n"
		<< L << "x_kp1 = x_k + d_k\n"
		<< L << "f_kp1 = f(x_kp1)\n"
		<< L << "c_kp1 = c(x_kp1)\n";
}
