// //////////////////////////////////////////////////////////////////
// print_algorithm_step.cpp
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

#include <ostream>

#include "IterationPack_print_algorithm_step.hpp"

void IterationPack::print_algorithm_step( const Algorithm& algo
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
