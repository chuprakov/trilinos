// ////////////////////////////////////////////////////////////////
// MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.h
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

#ifndef MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H

#include "MeritFunc_ModifiedL1LargerSteps_AddedStep.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{MeritFunc_ModifiedL1LargerSteps_AddedStep} from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group MeritFuncModifiedL1LargerSteps {
		after_k_iter                = 3;
		obj_increase_threshold      = 1e-3;
		max_pos_penalty_increase    = 1.0;
		pos_to_neg_penalty_increase = 1.0;
		incr_mult_factor            = 1e-4;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[after_k_iter] ToDo : Finish.
  *		Example: after_k_iter = 4;
  *	\item[obj_increase_threshold] ToDo : Finish.
  *		Example: obj_increase_threshold = 1e-4;
  *	\item[max_pos_penalty_increase] ToDo : Finish.
  *		Example: max_pos_penalty_increase = 1.0;
  *	\item[pos_to_neg_penalty_increase] ToDo : Finish.
  *		Example: pos_to_neg_penalty_increase = 1.0;
  *	\item[incr_mult_factor] ToDo : Finish.
  *		Example: incr_mult_factor = 1e-4;
  *	\end{description}
  */
class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			MeritFunc_ModifiedL1LargerSteps_AddedStep >
{
public:

	///
	MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions(
		MeritFunc_ModifiedL1LargerSteps_AddedStep* target = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
