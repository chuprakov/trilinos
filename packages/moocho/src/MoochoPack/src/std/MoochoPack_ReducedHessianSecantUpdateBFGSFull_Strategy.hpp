// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSFull_Strategy.hpp
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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H

#include "MoochoPack/src/std/ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack/src/std/BFGSUpdate_Strategy.hpp"
#include "MoochoPack/src/std/quasi_newton_stats.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Perform BFGS updates on full reduced Hessian.
 *
 * This is really a do nothing class that just uses a strategy
 * object (see #bfgs_update# below) to perform the update on the full
 * reduced hessian matrix.
 */
class ReducedHessianSecantUpdateBFGSFull_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
	
	///
	/** <<std comp>> members for the strategy object that will
	 * perform guts secant update.
	 */
	STANDARD_COMPOSITION_MEMBERS( BFGSUpdate_Strategy, bfgs_update )

    ReducedHessianSecantUpdateBFGSFull_Strategy(
		const bfgs_update_ptr_t&      bfgs_update = Teuchos::null
		);      

	/** @name Overridden from ReducedHessianSecantUpdate_Strategy */
	//@{
	///
	bool perform_update(
		VectorMutable     *s_bfgs
		,VectorMutable    *y_bfgs
		,bool                   first_update
		,std::ostream           & out
		,EJournalOutputLevel    olevel
		,NLPAlgo               *algo
		,NLPAlgoState              *s
		,MatrixSymOp        *rHL_k
		);
	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;
	//@}

private:
	quasi_newton_stats_iq_member	quasi_newton_stats_;

}; // end class ReducedHessianSecantUpdateBFGSFull_Strategy

}  // end namespace MoochoPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H
