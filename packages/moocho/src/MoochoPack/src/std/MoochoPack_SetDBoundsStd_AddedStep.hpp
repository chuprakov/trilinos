// ////////////////////////////////////////////////////////////////////////////
// SetDBoundsStd_AddedStep.h

#ifndef SET_D_BOUNDS_STD_ADDED_STEP_HH
#define SET_D_BOUNDS_STD_ADDED_STEP_HH

#include "ReducedSpaceSQPPack/include/rSQPAlgo_Step.h"
#include "ReducedSpaceSQPPack/include/std/d_bounds_iter_quant.h"

namespace ReducedSpaceSQPPack {

///
/** Updates the active set statistics for the current iteration.
  *
  * The function act_set_stats(...) is used to access the iteration
  * quantity object for the active set object.
  */
class SetDBoundsStd_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	d_bounds_iq_member		d_bounds_;

};	// end class SetDBoundsStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// SET_D_BOUNDS_STD_ADDED_STEP_HH