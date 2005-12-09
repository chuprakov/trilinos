// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_NewDecompositionSelectionStd_Strategy.hpp
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

#ifndef NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
#define NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Just force the decomposition system object to select a new
 * decomposition and let everyone else fend for themselves.
 */
class NewDecompositionSelectionStd_Strategy
	: public NewDecompositionSelection_Strategy
{
public:

	/// «std comp» members for range/null decomposition handler
	STANDARD_COMPOSITION_MEMBERS( DecompositionSystemHandlerSelectNew_Strategy, decomp_sys_handler )

	///
	NewDecompositionSelectionStd_Strategy(
		const decomp_sys_handler_ptr_t   &decomp_sys_handler
		);

	/** @name Overridden from NewDecompositionSelection_Strategy */
	//@{

	bool new_decomposition(
		NLPAlgo& algo, Algorithm::poss_type step_poss
		,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
		);
	///
	void print_new_decomposition(
		const NLPAlgo& algo, Algorithm::poss_type step_poss
		,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
		,std::ostream& out, const std::string& leading_str
		) const;

	//@}

private:

	// Not defined and not to be called
	NewDecompositionSelectionStd_Strategy();

};	// end class NewDecompositionSelectionStd_Strategy

}	// end namespace MoochoPack 

#endif	// NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
