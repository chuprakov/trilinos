// ////////////////////////////////////////////////////////////////////////////
// DecompositionSystemHandler_Strategy.hpp
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

#ifndef DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H
#define DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.hpp"
#include "GeneralIterationPack/src/Algorithm.hpp"

namespace ReducedSpaceSQPPack {

///
/** Interface for range/null decomposition handling.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemHandler_Strategy {
public:
	
	/** @name Public types */
	//@{

	///
	enum EDecompSysTesting { DST_DEFAULT, DST_TEST, DST_NO_TEST };
	///
	enum EDecompSysPrintLevel { DSPL_USE_GLOBAL, DSPL_LEAVE_DEFAULT };

	//@}

	///
	~DecompositionSystemHandler_Strategy() {}

	///
	/** Update the decomposition.
	 *
	 * This method may select a new decomposition (permuting the variables
	 * and constriants) and/or take control of the algorithm.
	 */
	virtual bool update_decomposition(
		rSQPAlgo                                &algo
		,rSQPState                              &s
		,NLPFirstOrderInfo                      &nlp
		,EDecompSysTesting                      decomp_sys_testing
		,EDecompSysPrintLevel                   decomp_sys_testing_print_level
		,bool                                   *new_decomp_selected
		) = 0;

	///
	/** Print the algorithm used for updating the decomposition.
	 */
	virtual void print_update_decomposition(
		const rSQPAlgo                          &algo
		,const rSQPState                        &s
		,std::ostream                           &out
		,const std::string                      &leading_spaces
		) const = 0;

}; // end class DecompositionSystemHandler_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H
