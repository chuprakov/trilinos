// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.h

#ifndef CHECK_CONVERGENCE_STD_ADDEDSTEP_H
#define CHECK_CONVERGENCE_STD_ADDEDSTEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Check for convergence.
  */
class CheckConvergenceStd_AddedStep : public rSQPAlgo_Step {
public:

	///
	CheckConvergenceStd_AddedStep(bool check_d = false)
		 : check_d_(check_d)
	{}

	/// Check ||d||inf < tolerance ?
	virtual void set_check_d(bool check_d)
	{	check_d_ = check_d; }

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	bool	check_d_;	// flag for if I should check ||d||inf
	
};	// end class CheckConvergenceStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// CHECK_CONVERGENCE_STD_ADDEDSTEP_H