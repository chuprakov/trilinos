// //////////////////////////////////////////////////////////////////////////
// MoochoPack_QuasiRangeSpaceStepStd_Strategy.hpp
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

#ifndef QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H

#include "MoochoPack_QuasiRangeSpaceStep_Strategy.hpp"

namespace MoochoPack {

///
/** Strategy class for computing a quasi-range-space step by solving
 * the approximate range space problem directly.
 */
class QuasiRangeSpaceStepStd_Strategy : public QuasiRangeSpaceStep_Strategy {
public:

	/** @name Overridden from QuasiRangeSpaceStep_Strategy */
	//@{

	///
	/** Solves the range space problem with the old decomposition at x_k.
	 *
	 * Solves:
	 \verbatim
	 vy = -inv(Gc_k'*Y_k) * c_xo
	 v = Y_k*vy
	 \endverbatim
	 */
 	 bool solve_quasi_range_space_step(
		std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
		,const Vector& xo, const Vector& c_xo, VectorMutable* v
	  	);

	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

	//@}

}; // end class QuasiRangeSpaceStepStd_Strategy

} // end namespace MoochoPack

#endif // QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H
