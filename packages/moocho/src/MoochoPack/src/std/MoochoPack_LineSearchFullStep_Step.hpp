// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStep_Step.h

#ifndef LINE_SEARCH_FULL_STEP_STEP_H
#define LINE_SEARCH_FULL_STEP_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptimizationPack/include/VariableBoundsTester.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Takes the full step x_kp1 = x_k + d_k (d_k = Ypy_k + Zpz_k).
  */
class LineSearchFullStep_Step : public LineSearch_Step {
public:

	///
	typedef ConstrainedOptimizationPack::VariableBoundsTester
		VariableBoundsTester;

	/// «std comp» Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )

	///
	LineSearchFullStep_Step(
		const bounds_tester_ptr_t&	bounds_tester
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
	/// Not defined and not to be called
	LineSearchFullStep_Step();

};	// end class LineSearchFullStep_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_FULL_STEP_STEP_H
