// ////////////////////////////////////////////////////////////////////////////////
// FeasibilityStep_Strategy.h
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

#ifndef FEASIBILITY_STEP_STRATEGY_H
#define FEASIBILITY_STEP_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Abstract interface for a strategy object that will compute a step that will
 * improve feasibility (at least descent) {abstract}.
 */
class FeasibilityStep_Strategy {
public:

	///
	virtual ~FeasibilityStep_Strategy() {}

	///
	/** Compute a step that improves feasibility (at least locally).
	 *
	 * This function will compute a step <tt>w</tt> that satisfies:
	 *
	 * <tt>d_bounds_k.l <= xo + w - x_k <= d_bounds_k.u</tt>
	 *
	 * and gives descent for <tt>||c(xo + beta*w)||</tt> for at least small <tt>beta > 0</tt>.
	 * This norm <tt>||.||</tt> could be any valid norm and the implementation is free to
	 * define what descent means any way it would like.  Any information being used
	 * in the algorithm can be used to compute this step.
	 *
	 * @param out     [out] Output stream journal data is written to.
	 * @param olevel  [in] Output level for printing to #out#
	 * @param algo    [in/out] The rSQPAlgo object.  This object can be queryed for
	 *                information.
	 * @param s       [in/out] rSQPState object.  May be queried or modified if needed.
	 * @param xo      [in] Base point vector (size n) xo.
	 * @param c_xo    [in] c(xo).
	 * @param w       [out] Computed step vector (size n) w.  Must not be NULL.
	 *
	 * @return Returns true if a step that reduces feasibility subject to the bounds could
	 * be found and false otherwise.
	 */
 	virtual bool compute_feasibility_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorWithOp& xo, const VectorWithOp& c_xo, VectorWithOpMutable* w
	  	) = 0;

	///
	/** This function will print a description of the computations and logic used.
	 */
	virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class FeasibilityStep_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // FEASIBILITY_STEP_STRATEGY_H
