// ////////////////////////////////////////////////////////////////
// NLPInterfacePack_NLPFirstDerivTesterSetOptions.hpp
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

#ifndef NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H
#define NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H

#include "NLPInterfacePack_NLPFirstDerivTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

///
/** Set options for NLPFirstDerivTester from an
 * OptionsFromStream object.
 *
 * The default options group name is NLPFirstDerivTester.
 *
 * The options group is:
 *
 \verbatim
 options_group NLPFirstDerivTester {
 *    fd_testing_method = FD_COMPUTE_ALL;
     fd_testing_method = FD_DIRECTIONAL;
     num_fd_directions = 3;  *** [fd_testing_method == DIRECTIONAL]
     warning_tol   = 1e-14;
     error_tol     = 1e-8;
 }
 \endverbatim
 */
class NLPFirstDerivTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPFirstDerivTester >
{
public:

	///
	NLPFirstDerivTesterSetOptions(
		  NLPFirstDerivTester* target = 0
		, const char opt_grp_name[] = "NLPFirstDerivTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class NLPFirstDerivTesterSetOptions

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H
