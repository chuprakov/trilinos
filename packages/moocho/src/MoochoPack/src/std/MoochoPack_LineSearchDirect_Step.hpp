// ////////////////////////////////////////////////////////////////////////////
// LineSearchDirect_Step.h
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

#ifndef LINE_SEARCH_DIRECT_STEP_H
#define LINE_SEARCH_DIRECT_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearch_Strategy.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardAggregationMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Delegates the line search (d_k = Ypy_k + Zpz_k).
  */
class LineSearchDirect_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for direct_line_search
	STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search)

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	LineSearchDirect_Step()
		: direct_line_search_(0), merit_func_(0)
	{}

	///
	LineSearchDirect_Step(const direct_line_search_ptr_t& direct_line_search
			, const merit_func_ptr_t& merit_func)
		: direct_line_search_(direct_line_search), merit_func_(merit_func)
	{}

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class LineSearchDirect_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_DIRECT_STEP_H
