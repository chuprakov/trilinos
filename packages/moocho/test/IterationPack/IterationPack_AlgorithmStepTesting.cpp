// ///////////////////////////////////////////////////////////////
// AlgorithmStepTesting.cpp
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

#include <iomanip>

#include "AlgorithmStepTesting.h"
#include "GeneralIterationPack/src/Algorithm.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"

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
