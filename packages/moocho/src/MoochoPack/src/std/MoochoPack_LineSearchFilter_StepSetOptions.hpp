// ////////////////////////////////////////////////////////////////
// MoochoPack_LineSearchFilter_StepSetOptions.hpp
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

#ifndef LINE_SEARCH_FILTER_STEP_SET_OPTIONS_H
#define LINE_SEARCH_FILTER_STEP_SET_OPTIONS_H

#include "MoochoPack_LineSearchFilter_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for LineSearchFilter_Step from an \c OptionsFromStream object.
 *
 * The default options group name is LineSearchFilter.
 *
 * The options group is:
 *
 \verbatim
 options_group LineSearchFilter {
 gamma_theta     = 0.01;
 gamma_psi       = 0.01;
 back_track_frac = 0.5;
 }
 \verbatim
*/
class LineSearchFilter_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode,
	  public OptionsFromStreamPack::SetOptionsToTargetBase< LineSearchFilter_Step >
	{
	public:

		///
		LineSearchFilter_StepSetOptions(
		  LineSearchFilter_Step* target = 0,
		  const char opt_grp_name[] = "LineSearchFilter" );

	protected:

		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );
	
	};	// end class LineSearchFilter_StepSetOptions

}	// end namespace MoochoPack

#endif	// LINE_SEARCH_FILTER_STEP_SET_OPTIONS_H
