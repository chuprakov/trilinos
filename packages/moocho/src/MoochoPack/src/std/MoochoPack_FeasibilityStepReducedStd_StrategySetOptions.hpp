// ////////////////////////////////////////////////////////////////
// FeasibilityStepReducedStd_StrategySetOptions.h

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H

#include "FeasibilityStepReducedStd_Strategy.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for FeasibilityStepReducedStd_Strategy from an
  * OptionsFromStream object.
  *
  * The default options group name is IndepDirecWithBoundsStd.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group FeasibilityStepReducedStd_Strategy {
	*    qp_objective = OBJ_MIN_FULL_STEP;
	*    qp_objective = OBJ_MIN_NULL_SPACE_STEP;
	*    qp_objective = OBJ_RSQP;
	*    qp_testing   = QP_TEST_DEFAULT;
	*    qp_testing   = QP_TEST;
	*    qp_testing   = QP_NO_TEST;
	}
  \end{verbatim}
  */
class FeasibilityStepReducedStd_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			FeasibilityStepReducedStd_Strategy >
{
public:

	///
	FeasibilityStepReducedStd_StrategySetOptions(
		  FeasibilityStepReducedStd_Strategy* target = 0
		, const char opt_grp_name[] = "FeasibilityStepReducedStd" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class FeasibilityStepReducedStd_StrategySetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
