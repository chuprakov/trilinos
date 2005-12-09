// ////////////////////////////////////////////////////////////////
// NLPInterfacePack_NLPDirectTesterSetOptions.hpp
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

#ifndef NLP_FIRST_ORDER_DIRECT_TESTER_SET_OPTIONS_H
#define NLP_FIRST_ORDER_DIRECT_TESTER_SET_OPTIONS_H

#include "NLPInterfacePack_NLPDirectTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

///
/** Set options for <tt>NLPDirectTester</tt> from an
  * <tt>\ref OptionsFromStreamPack::OptionsFromStream "OptionsFromStream"</tt> object.
  *
  * The default options group name is 'NLPDirectTester'.
  *
  * The options group is:
  \verbatim

	options_group NLPDirectTester {
	*    Gf_testing_method = FD_COMPUTE_ALL;
	    Gf_testing_method = FD_DIRECTIONAL;
	    Gf_warning_tol    = 1e-6;
	    Gf_error_tol      = 1e-1;
	*    Gc_testing_method = FD_COMPUTE_ALL;
	    Gc_testing_method = FD_DIRECTIONAL;
	    Gc_warning_tol    = 1e-6;
	    Gc_error_tol      = 1e-1;
		num_fd_directions = 3;  *** [testing_method == FD_DIRECTIONAL]
	}
  \endverbatim
  */
class NLPDirectTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPDirectTester >
{
public:

	///
	NLPDirectTesterSetOptions(
		  NLPDirectTester* target = 0
		, const char opt_grp_name[] = "NLPDirectTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class NLPDirectTesterSetOptions

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_ORDER_DIRECT_TESTER_SET_OPTIONS_H
