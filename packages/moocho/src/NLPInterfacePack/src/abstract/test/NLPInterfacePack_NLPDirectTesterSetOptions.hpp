// ////////////////////////////////////////////////////////////////
// NLPFirstOrderDirectTesterSetOptions.hpp
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

#include "NLPFirstOrderDirectTester.hpp"
#include "SetOptionsFromStreamNode.hpp"
#include "SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

///
/** Set options for <tt>NLPFirstOrderDirectTester</tt> from an
  * <tt>\ref OptionsFromStreamPack::OptionsFromStream "OptionsFromStream"</tt> object.
  *
  * The default options group name is 'NLPFirstOrderDirectTester'.
  *
  * The options group is:
  \verbatim

	options_group NLPFirstOrderDirectTester {
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
class NLPFirstOrderDirectTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPFirstOrderDirectTester >
{
public:

	///
	NLPFirstOrderDirectTesterSetOptions(
		  NLPFirstOrderDirectTester* target = 0
		, const char opt_grp_name[] = "NLPFirstOrderDirectTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class NLPFirstOrderDirectTesterSetOptions

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_ORDER_DIRECT_TESTER_SET_OPTIONS_H
