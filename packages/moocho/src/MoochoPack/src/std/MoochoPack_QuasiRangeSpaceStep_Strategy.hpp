// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStep_Strategy.h
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

#ifndef QUASI_RANGE_SPACE_STEP_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Abstract interface for a strategy object that will compute a step that will
 * approximalty solve a range space subproblem {abstract}.
 */
class QuasiRangeSpaceStep_Strategy {
public:

	///
	virtual ~QuasiRangeSpaceStep_Strategy() {}

	///
	/** Compute a step that will approximatly solve a range space subproblem.
	 *
	 * This function will compute a step #v# that will approximatly satisfy:
	 *
	 * #||Gc_k'*v + c(xo)|| < ||c(xo)||#\\
	 *
	 * The above norm ||.|| could be any valid norm and the implementation is free to
	 * define what descent means any way it would like.  It is assumed that this
	 * step will be computed by using the #Gc_k# but other implementations are
	 * possible.  Any information being used in the algorithm can be used
	 * to compute this step in a reasonable way.  Note that the 
	 * inequalities do not have to (and should not in most cases) be considered
	 * in this computation.  Note that whatever means is used to compute #v#
	 * that it better give a descent direction for #||c(x)||# but there is no guarantee
	 * for this if #||xo - x_k||# is large since #Gc_k# may not accurately approximate
	 * #Gc(xo)#.
	 *
	 * @param out      [out] Output stream journal data is written to.
	 * @param olevel   [in] Output level for printing to #out#.
	 * @param algo     [in/out] The rSQPAlgo object.  This object can be queryed for information.
	 * @param s        [in/out] rSQPState object.  May be queried or modified if needed.
	 * @param xo       [in] Base point vector (size n) xo.
	 * @param c_xo     [out] Constraints residual c(xo).
	 * @param v        [out] Computed step vector (size n).  Must not be NULL.
	 *
	 * @return Returns true if a step could be found and false otherwise.
	 */
 	virtual bool solve_quasi_range_space_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* v
	  	) = 0;

	///
	/** This function will print a description of the computations and logic used.
	 */
	virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class QuasiRangeSpaceStep_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // QUASI_RANGE_SPACE_STEP_STRATEGY_H
