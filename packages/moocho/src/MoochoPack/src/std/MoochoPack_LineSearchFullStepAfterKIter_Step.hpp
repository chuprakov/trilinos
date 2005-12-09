// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_LineSearchFullStepAfterKIter_Step.hpp
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

#ifndef LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
#define LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H

#include <limits>

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"

namespace MoochoPack {

///
/** Changes from a line search step to just taking full steps after
  * full_steps_after_k iterations.
  */
class LineSearchFullStepAfterKIter_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for the line search step
	STANDARD_COMPOSITION_MEMBERS(LineSearch_Step,line_search)

	///
	LineSearchFullStepAfterKIter_Step(
			  const line_search_ptr_t&	line_search			= 0
			, int						full_steps_after_k
											= std::numeric_limits<int>::max()	)
		: line_search_(line_search)
			, full_steps_after_k_(full_steps_after_k)
	{}

	/// 
	void full_steps_after_k( int full_steps_after_k )
	{	full_steps_after_k_ = full_steps_after_k; }
	///
	value_type full_steps_after_k() const
	{	return full_steps_after_k_; }

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	int		full_steps_after_k_;
	
};	// end class LineSearchFullStepAfterKIter_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
