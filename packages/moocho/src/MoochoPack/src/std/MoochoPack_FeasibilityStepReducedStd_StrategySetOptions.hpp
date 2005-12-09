// ////////////////////////////////////////////////////////////////
// MoochoPack_FeasibilityStepReducedStd_StrategySetOptions.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_FeasibilityStepReducedStd_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

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
	void setOption( int option_num, const std::string& option_value );

};	// end class FeasibilityStepReducedStd_StrategySetOptions

}	// end namespace MoochoPack

#endif	// FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
