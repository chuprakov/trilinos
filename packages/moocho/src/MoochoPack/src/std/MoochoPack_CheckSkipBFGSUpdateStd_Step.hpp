// ////////////////////////////////////////////////////////////////////////////
// CheckSkipBFGSUpdateStd_Step.h

#ifndef CHECK_SKIP_BFGS_UPDATE_STD_STEP_H
#define CHECK_SKIP_BFGS_UPDATE_STD_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"

namespace ReducedSpaceSQPPack {

///
/** Checks if a BFGS update should be preformed.
  */
class CheckSkipBFGSUpdateStd_Step : public rSQPAlgo_Step {
public:

	///
	CheckSkipBFGSUpdateStd_Step()
	{}

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	
};	// end class ReducedHessianBFGS_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// CHECK_SKIP_BFGS_UPDATE_STD_STEP_H