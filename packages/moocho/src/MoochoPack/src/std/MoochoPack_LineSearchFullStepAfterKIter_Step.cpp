// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStepAfterKIter_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/LineSearchFullStepAfterKIter_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"

bool ReducedSpaceSQPPack::LineSearchFullStepAfterKIter_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLP			&nlp	= algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	out << std::boolalpha;

	// print step header.
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	const bool take_full_step = s.k() > full_steps_after_k();

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "\nk = " << s.k() << ( take_full_step ? " > " : " < ")
					<< "full_steps_after_k = " << full_steps_after_k() << std::endl;
	}

	if( !take_full_step ) {
		return line_search().do_step(_algo,step_poss,type,assoc_step_poss);
	}
	else {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "\nKeep the full step...\n";
		}
	}

	return true;
}

void ReducedSpaceSQPPack::LineSearchFullStepAfterKIter_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out	<< L << "*** Start using full steps after full_steps_after_k iterations.\n"
		<< L << "default: full_steps_after_k = very big\n";
	out	<< L << "if k < full_steps_after_k then\n";
	line_search().print_step(algo,step_poss,type,assoc_step_poss,out,L + "    " );
	out	<< L << "end\n";
}
