// ///////////////////////////////////////////////////////
// MoochoPack_ReducedHessianSecantUpdate_Strategy.hpp
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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

///
/** Strategy interface for performing secant updates {abstract}.
 *
 * This interface is used by the class \c ReducedHessianSecantUpdateStd_Step
 * to actually perform the secant updates.
 */
class ReducedHessianSecantUpdate_Strategy {
public:
	
    ///
	virtual ~ReducedHessianSecantUpdate_Strategy() {}

	///
	/** Perform the secant update.
	 *
	 * The function will update <tt>rHL_k</tt> so that <tt>rHL_k * s_bfgs \approx y_bfgs</tt>.
	 * Note that this post conditions for this function do not strictly require
	 * that the secant property <tt>rHL_k * s_bfgs = y_bfgs</tt> be satisfied.  This
	 * allows for more flexibility in how the update is perform.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>s_bfgs->size() == y_bfgs->size() == rHL_k->rows() == rHL_k->cols()</tt> (throws ???)
	 * </ul>
	 *
	 * @param  s_bfgs [in/work] Secant change vector on input.  May be modified as
	 *                modified as workspace.
	 * @param  y_bfgs [in/work] Secant change vector on input.  May be modified as
 	 *                modified as workspace.
	 * @param  first_update
	 *                [in] If true then this is the first update after <tt>rHL</tt> was
	 *                initialized to identity.  This is information that may be
	 *                used in order to deliver a beter initial update.
	 * @param  out    [out] Output stream journal data is written to.
	 * @param  olevel [in] Output level for printing to <tt>out</tt>
	 * @param  algo   [in/out] The NLPAlgo object.  This object can be queryed for
	 *                information and also be called to redirect control (in which
	 *                case this function should probably return false).
	 * @param  s      [in/out] NLPAlgoState object.  May be queried or modified if needed.
	 * @param  rHL_k  [in/out] The matrix to be updated.  Note that <tt>rHL_k</tt> was already
	 *                set to <tt>rHL_km1</tt> before this call was made.  Also, <tt>rHL_k</tt> will
	 *                probably have to support the <tt>MatrixSymSecant</tt> interface
	 *                or an exception will be thrown.
	 * 
	 * @return Returns false if the algorithms path has been redirected through <tt>algo</tt>.
	 * Ohterwise, this function should return true.
	 */
	virtual bool perform_update(
		VectorMutable     *s_bfgs
		,VectorMutable   *y_bfgs
		,bool                   first_update
		,std::ostream           & out
		,EJournalOutputLevel    olevel
		,NLPAlgo               *algo
		,NLPAlgoState              *s
		,MatrixSymOp        *rHL_k
		) = 0;
	
	///
	/** This function will print a description of the computations and logic used
	 * in the update.
	 */
	virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class ReducedHessianSecantUpdate_Strategy

}  // end namespace MoochoPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H
