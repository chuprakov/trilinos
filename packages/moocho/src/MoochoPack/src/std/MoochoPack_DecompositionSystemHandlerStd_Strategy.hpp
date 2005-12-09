// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_DecompositionSystemHandlerStd_Strategy.hpp
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

#ifndef DECOMPOSITION_SYSTEM_HANDLER_STD_STRATEGY_H
#define DECOMPOSITION_SYSTEM_HANDLER_STD_STRATEGY_H

#include "MoochoPack_DecompositionSystemHandler_Strategy.hpp"

namespace MoochoPack {

///
/** Subclass for updating the range/null space decomposition using
 * the base DecompositionSystem interface only.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemHandlerStd_Strategy
	: public DecompositionSystemHandler_Strategy
{
public:
	
	/** @name Constructors / initializers */
	//@{

	///
	/** Constructor
	 */
	DecompositionSystemHandlerStd_Strategy();

	//@}

	/** @name Overridden from DecompositionSystemHandler_Strategy */
	//@{

	///
	bool update_decomposition(
		NLPAlgo                                &algo
		,NLPAlgoState                          &s
		,NLPFirstOrder                         &nlp
		,EDecompSysTesting                     decomp_sys_testing
		,EDecompSysPrintLevel                  decomp_sys_testing_print_level
		,bool                                  *new_decomp_selected
		);
	///
	void print_update_decomposition(
		const NLPAlgo                          &algo
		,const NLPAlgoState                    &s
		,std::ostream                          &out
		,const std::string                     &leading_spaces
		) const;

	//@}

}; // end class DecompositionSystemHandlerStd_Strategy

} // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_HANDLER_STD_STRATEGY_H
