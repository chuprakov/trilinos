// ////////////////////////////////////////////////////////////////
// MoochoPack_CheckSkipBFGSUpdateStd_StepSetOptions.hpp
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

#ifndef CHECK_SKIP_BFGS_UPDATE_STD_STEP_SET_OPTIONS_H
#define CHECK_SKIP_BFGS_UPDATE_STD_STEP_SET_OPTIONS_H

#include "MoochoPack_CheckSkipBFGSUpdateStd_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for CheckSkipBFGSUpdateStd_Step from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group CheckSkipBFGSUpdateStd {
		skip_bfgs_prop_const	= 10.0; *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[skip_bfgs_prop_const] ToDo : Finish.
  *		Example: skip_bfgs_prop_const = 10.0;
  *	\end{description}
  */
class CheckSkipBFGSUpdateStd_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			CheckSkipBFGSUpdateStd_Step >
{
public:

	///
	CheckSkipBFGSUpdateStd_StepSetOptions(
		CheckSkipBFGSUpdateStd_Step* target = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class CheckSkipBFGSUpdateStd_StepSetOptions

}	// end namespace MoochoPack

#endif	// CHECK_SKIP_BFGS_UPDATE_STD_STEP_SET_OPTIONS_H
