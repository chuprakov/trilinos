// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepStd_Strategy.h

#ifndef QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H

#include "QuasiRangeSpaceStep_Strategy.h"

namespace ReducedSpaceSQPPack {

///
/** Strategy class for computing a quasi range space step by solving
 * the approximate range space problem directly.
 */
class QuasiRangeSpaceStepStd_Strategy : public QuasiRangeSpaceStep_Strategy {
public:

	// ////////////////////////////////////////////////////////////
	// Overridden from QuasiRangeSpaceStep_Strategy

	///
	/** Solves the range space problem with the old decomposition at x_k.
	 *
	 * Solves:
	 *
	 * #vy = -inv(Gc_k'*Y_k) * c_xo#
	 * #v = Y_k*vy#
	 */
 	 bool solve_quasi_range_space_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* v
	  	);

	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

}; // end class QuasiRangeSpaceStepStd_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // QUASI_RANGE_SPACE_STEP_STD_STRATEGY_H
