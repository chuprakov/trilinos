// ////////////////////////////////////////////////////////////////////////////
// CrossTermExactStd_Step.h

#ifndef CROSS_TERM_EXACT_STD_STEP_H
#define CROSS_TERM_EXACT_STD_STEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** w_k = Z_k' * HL_k * Ypy_k
  */
class CrossTermExactStd_Step : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CrossTermExactStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// CROSS_TERM_EXACT_STD_STEP_H
