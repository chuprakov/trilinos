// ////////////////////////////////////////////////////////////////
// ConstrainedOptPack_VariableBoundsTesterSetOptions.hpp
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

#ifndef VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H
#define VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H

#include "ConstrainedOptPack_VariableBoundsTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

///
/** Set options for VariableBoundsTester from an OptionsFromStream
  * object.
  *
  * The default options group name is VariableBoundsTester.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group VariableBoundsTester {
	    warning_tol   = 1e-10;
	    error_tol     = 1e-5;
	}
  \end{verbatim}
  */
class VariableBoundsTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			VariableBoundsTester >
{
public:

	///
	VariableBoundsTesterSetOptions(
		  VariableBoundsTester* target = 0
		, const char opt_grp_name[] = "VariableBoundsTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class VariableBoundsTesterSetOptions

}	// end namespace ConstrainedOptPack

#endif	// VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H
