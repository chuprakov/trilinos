// ///////////////////////////////////////////////////////////////
// AlgorithmStepTesting.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <iomanip>

#include "../test/AlgorithmStepTesting.h"
#include "../include/Algorithm.h"
#include "../include/print_algorithm_step.h"

namespace GeneralIterationPack {

bool AlgorithmStepTesting::do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss)
{
	print_algorithm_step( algo, step_poss, type, assoc_step_poss
		, algo.track().journal_out() );
	return true;
}

void AlgorithmStepTesting::inform_updated(Algorithm& algo)
{
	algo.track().journal_out() << "\ninform_updated(algo) called for step\n";
}

void AlgorithmStepTesting::print_step( const Algorithm& algo, poss_type step_poss, EDoStepType type
	, poss_type assoc_step_poss ,std::ostream& out, const std::string& leading_str ) const
{
	char type_name[3][15] = { "DO_MAIN_STEP", "DO_PRE_STEP" , "DO_POST_STEP" };
	algo.track().journal_out()
		<< std::endl << leading_str << step_poss << ", " << type_name[type];
	if(type == DO_MAIN_STEP) {
		algo.track().journal_out()
			<< ", \"" << algo.get_step_name(step_poss) << "\"";
	}
	else {
		EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
		algo.track().journal_out()
			<< ", " << assoc_step_poss
			<< ", \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
	}
	algo.track().journal_out()
		<< "\" : print_step(algo,step_poss,type,assoc_step_poss,out) called\n";
}

}	// end namespace GeneralIterationPack 
