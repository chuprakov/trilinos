// ////////////////////////////////////////////////////////////////////////////
// CheckSkipBFGSUpdateStd_Step.h

#ifndef CHECK_SKIP_BFGS_UPDATE_STD_STEP_H
#define CHECK_SKIP_BFGS_UPDATE_STD_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Checks if a BFGS update should be preformed.
  */
class CheckSkipBFGSUpdateStd_Step : public rSQPAlgo_Step {
public:

	///
	/** <<std member comp>> members for proportionality constant to use in the
	  * test for if to perform BFGS update.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, skip_bfgs_prop_const )

	///
	CheckSkipBFGSUpdateStd_Step(
		value_type	skip_bfgs_prop_const	= 10.0
		);

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