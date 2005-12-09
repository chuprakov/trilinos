// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp
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

#ifndef DECOMPOSITION_SYSTEM_HANDLER_SELECT_NEW_STRATEGY_H
#define DECOMPOSITION_SYSTEM_HANDLER_SELECT_NEW_STRATEGY_H

#include "MoochoPack_DecompositionSystemHandler_Strategy.hpp"

namespace MoochoPack {

///
/** Interface for range/null decomposition handling.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemHandlerSelectNew_Strategy
	: public DecompositionSystemHandler_Strategy
{
public:
	
	///
	/** Instruct the \c DecompositionSystemHandler_Strategy object to select a new decomposition
	 * the next time \c update_decomposition() is called.
	 */
	virtual void select_new_decomposition( bool select_new_decomposition = true ) = 0;

}; // end class DecompositionSystemHandlerSelectNew_Strategy

} // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_HANDLER_SELECT_NEW_STRATEGY_H
