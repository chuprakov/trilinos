// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepTailoredApproach_Strategy.h

#ifndef QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H

#include "QuasiRangeSpaceStep_Strategy.h"

namespace ReducedSpaceSQPPack {

///
/** Strategy class for computing a quasi range space step for the
 * tailored approach NLP interface.
 */
class QuasiRangeSpaceStepTailoredApproach_Strategy : public QuasiRangeSpaceStep_Strategy {
public:

	// ////////////////////////////////////////////////////////////
	// Overridden from QuasiRangeSpaceStep_Strategy

	///
	/** Calls the NLPrSQPTailoredApproach iterface to compute the step.
	 *
	 * ToDo: Finish documentation!
	 */
 	 bool solve_quasi_range_space_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* v
	  	);

	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

}; // end class QuasiRangeSpaceStepTailoredApproach_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H
