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

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/src/AlgorithmStep.h"
#include "ConstrainedOptimizationPack/src/DirectLineSearch_Strategy.h"
#include "StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Delegates the line search to a <tt>DirectLineSearch_Strategy</tt> object.
 */
class LineSearchDirect_Step
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/// <<std comp>> members for direct_line_search
	STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search)
	///
	LineSearchDirect_Step(
		const direct_line_search_ptr_t& direct_line_search = MemMngPack::null
		);

	/** Overridden from AlgorithmStep */
	//@{

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

	//@}

};	// end class LineSearchDirect_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_DIRECT_STEP_H
