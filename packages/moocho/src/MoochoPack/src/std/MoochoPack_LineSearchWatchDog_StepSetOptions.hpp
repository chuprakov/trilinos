// ////////////////////////////////////////////////////////////////
// LineSearchWatchDog_StepSetOptions.hpp
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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
#define LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H

#include "LineSearchWatchDog_Step.hpp"
#include "SetOptionsFromStreamNode.hpp"
#include "SetOptionsToTargetBase.hpp"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{LineSearchWatchDog_Step} from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group LineSearchWatchDog {
		opt_kkt_err_threshold	= 1e-3; *** (+dbl)
		feas_kkt_err_threshold	= 1e-3; *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[opt_kkt_err_threshold] ToDo : Finish.
  *		Example: opt_kkt_err_threshold = 1e-1;
  *	\item[feas_kkt_err_threshold] ToDo : Finish.
  *		Example: feas_kkt_err_threshold = 1e-2;
  *	\end{description}
  */
class LineSearchWatchDog_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			LineSearchWatchDog_Step >
{
public:

	///
	LineSearchWatchDog_StepSetOptions(
		LineSearchWatchDog_Step* target = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class LineSearchWatchDog_StepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
