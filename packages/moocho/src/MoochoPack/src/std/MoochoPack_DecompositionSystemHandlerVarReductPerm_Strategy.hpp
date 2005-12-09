// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_DecompositionSystemHandlerVarReductPerm_Strategy.hpp
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

#ifndef DECOMPOSITION_SYSTEM_HANDLER_VAR_REDUCT_PERM_STRATEGY_H
#define DECOMPOSITION_SYSTEM_HANDLER_VAR_REDUCT_PERM_STRATEGY_H

#include "MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp"

namespace MoochoPack {

///
/** Subclass for selecting and updating the range/null space decomposition using
 * the DecompositionSystemVarReductPerm interface.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemHandlerVarReductPerm_Strategy
	: public DecompositionSystemHandlerSelectNew_Strategy
{
public:
	
	/** @name Constructors / initializers */
	//@{

	///
	/** Constructor
	 */
	DecompositionSystemHandlerVarReductPerm_Strategy();

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

	/** @name Overridden from DecompositionSystemHandlerSelectNew_Strategy */
	//@{

	///
	void select_new_decomposition( bool select_new_decomposition );

	//@}

private:

	bool select_new_decomposition_;

}; // end class DecompositionSystemHandlerVarReductPerm_Strategy

} // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_HANDLER_VAR_REDUCT_PERM_STRATEGY_H
