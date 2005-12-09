// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_ActSetStats_AddedStep.hpp
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

#ifndef ACT_SET_STATS_ADDED_STEP_H
#define ACT_SET_STATS_ADDED_STEP_H

#include "MoochoPack/src/rSQPAlgo_Step.h"
#include "MoochoPack_act_set_stats.hpp"

namespace MoochoPack {

///
/** Updates the active set statistics for the current iteration.
  *
  * The function act_set_stats(...) is used to access the iteration
  * quantity object for the active set object.
  */
class ActSetStats_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	act_set_stats_iq_member		act_set_stats_;

};	// end class ActSetStats_AddedStep

}	// end namespace MoochoPack 

#endif	// ACT_SET_STATS_ADDED_STEP_H
