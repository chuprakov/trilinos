// ////////////////////////////////////////////////////////////////
// MoochoPack_EvalNewPointTailoredApproach_StepSetOptions.hpp
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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for EvalNewPointTailoredApproach_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is EvalNewPointTailoredApproach.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group EvalNewPointTailoredApproach {
		fd_deriv_testing   = FD_DEFAULT;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[fd_deriv_testing] Determines if finite differerece testing of the 
  *		derivatives of the Gc and Gf.  See the class \Ref{EvalNewPointTailoredApproach_Step}
  *		and its printed algorithm for more details).
  *		\begin{description}
  *		\item[FD_DEFAULT]			The global flag check_results determines
  *									if the tests are performed.
  *		\item[FD_TEST]				The tests are performed reguardless the
  *									value of check_results
  *		\item[FD_NO_TEST]			The tests are not performed reguardless the
  *									value of check_results
  *		\end{description}
  */
class EvalNewPointTailoredApproach_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			EvalNewPointTailoredApproach_Step >
{
public:

	///
	EvalNewPointTailoredApproach_StepSetOptions(
		  EvalNewPointTailoredApproach_Step* target = 0
		, const char opt_grp_name[] = "EvalNewPointTailoredApproach" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class EvalNewPointTailoredApproach_StepSetOptions

}	// end namespace MoochoPack

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
