// ////////////////////////////////////////////////////////////////////////////
// SetDBoundsStd_AddedStep.h
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

#ifndef SET_D_BOUNDS_STD_ADDED_STEP_HH
#define SET_D_BOUNDS_STD_ADDED_STEP_HH

#include "ReducedSpaceSQPPack/include/rSQPAlgo_Step.h"
#include "ReducedSpaceSQPPack/include/std/d_bounds_iter_quant.h"

namespace ReducedSpaceSQPPack {

///
/** Updates the active set statistics for the current iteration.
  *
  * The function act_set_stats(...) is used to access the iteration
  * quantity object for the active set object.
  */
class SetDBoundsStd_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	d_bounds_iq_member		d_bounds_;

};	// end class SetDBoundsStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// SET_D_BOUNDS_STD_ADDED_STEP_HH
