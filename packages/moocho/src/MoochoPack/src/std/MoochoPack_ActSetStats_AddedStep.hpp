// ////////////////////////////////////////////////////////////////////////////
// ActSetStats_AddedStep.h

#ifndef ACT_SET_STATS_ADDED_STEP_H
#define ACT_SET_STATS_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Updates the active set statistics for the current iteration.
  *
  * The function act_set_stats(...) is used to access the iteration
  * quantity object for the active set object.
  */
class ActSetStats_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class ActSetStats_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// ACT_SET_STATS_ADDED_STEP_H