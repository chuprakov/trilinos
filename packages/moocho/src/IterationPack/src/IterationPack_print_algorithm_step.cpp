// //////////////////////////////////////////////////////////////////
// print_algorithm_step.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <ostream>

#include "../include/print_algorithm_step.h"

void GeneralIterationPack::print_algorithm_step( const Algorithm& algo
	, Algorithm::poss_type step_poss, EDoStepType type
	, Algorithm::poss_type assoc_step_poss, std::ostream& out )
{
	out << "\n(" << algo.state().k() << ") " << step_poss;
	if(type == DO_MAIN_STEP) {
		out
			<< ": \"" << algo.get_step_name(step_poss) << "\"";
	}
	else {
		EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
		int num_assoc_steps = algo.num_assoc_steps(step_poss,_type);
		out << ".";
		switch(_type) {
			case PRE_STEP:
				out << - num_assoc_steps + ((int)assoc_step_poss - 1);
				break;
			case POST_STEP:
				out << assoc_step_poss;
				break;
		}
		out	<< ": \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
	}
	out << "\n";
}