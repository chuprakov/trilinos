// ////////////////////////////////////////////////////////////////
// NullSpaceStepWithInequStd_StepSetOptions.h
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

#ifndef INDEP_DIREC_WITH_BOUNDS_STD_STEP_SET_OPTIONS_H
#define INDEP_DIREC_WITH_BOUNDS_STD_STEP_SET_OPTIONS_H

#include "NullSpaceStepWithInequStd_Step.h"
#include "SetOptionsFromStreamNode.h"
#include "SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for NullSpaceStepWithInequStd_Step from an
 * OptionsFromStream object.
 *
 * The default options group name is IndepDirecWithBoundsStd.
 *
 * The options group is:
 *
 \verbatim

	options_group NullSpaceStepWithInequStd {
	*    warm_start_frac   = 0.8;    *** (+dbl, [0.0,1.0]) Permform warm start when warm_start_frac * 100%
	                                 *** of the active set is the same.
	    qp_testing   = QP_TEST_DEFAULT;
	*    qp_testing   = QP_TEST;
	*    qp_testing   = QP_NO_TEST;
	*    primal_feasible_point_error = true;
	*    dual_feasible_point_error   = true;
	}
 \endverbatim
 */
class NullSpaceStepWithInequStd_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
	, public OptionsFromStreamPack::SetOptionsToTargetBase<
		NullSpaceStepWithInequStd_Step >
{
public:
	///
	NullSpaceStepWithInequStd_StepSetOptions(
		 NullSpaceStepWithInequStd_Step* target = NULL
		,const char opt_grp_name[] = "NullSpaceStepWithInequStd" );
protected:
	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );
};	// end class NullSpaceStepWithInequStd_StepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// INDEP_DIREC_WITH_BOUNDS_STD_STEP_SET_OPTIONS_H
