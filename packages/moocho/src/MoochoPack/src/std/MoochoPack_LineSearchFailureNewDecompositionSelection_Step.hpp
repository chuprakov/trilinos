// ////////////////////////////////////////////////////////////////////////////
// LineSearchFailureNewDecompositionSelection_Step.hpp
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

#ifndef LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
#define LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H

#include "NewDecompositionSelection_Strategy.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Directs the selection of a new decomposition if the line search fails.
  *
  * If the delegated line search Step object throws a \c LineSearchFailure
  * exception, then this object directs the selection of a new
  * decomposition.  If the very next iteration also results in a linesearch
  * failure then we must quit.
  */
class LineSearchFailureNewDecompositionSelection_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/// <<std comp>> members for LineSearch object.
	STANDARD_COMPOSITION_MEMBERS( IterationPack::AlgorithmStep, line_search_step )

	/// <<std comp>> members for Decomposition Select Strategy object.
	STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy )

	///
	LineSearchFailureNewDecompositionSelection_Step(
		const line_search_step_ptr_t        &line_search_step
		,const new_decomp_strategy_ptr_t    &new_decomp_strategy
		);

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:

	int last_ls_failure_k_;

	// not defined and not to be called
	LineSearchFailureNewDecompositionSelection_Step();
	LineSearchFailureNewDecompositionSelection_Step(
		const LineSearchFailureNewDecompositionSelection_Step&);
	LineSearchFailureNewDecompositionSelection_Step& operator=(
		const LineSearchFailureNewDecompositionSelection_Step&);

};	// end class LineSearchFailureNewDecompositionSelection_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
