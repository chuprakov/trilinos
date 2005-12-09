// ////////////////////////////////////////////////////////////////
// NLPInterfacePack_NLPTesterSetOptions.hpp
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

#ifndef NLP_TESTER_SET_OPTIONS_H
#define NLP_TESTER_SET_OPTIONS_H

#include "NLPInterfacePack_NLPTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

///
/** Set options for <tt>NLPTester</tt> from an
 * <tt>\ref OptionsFromStreamPack::OptionsFromStream "OptionsFromStream"</tt> object.
 *
 * The default options group name is '%NLPTester'.
 *
 * The options group is:
 \verbatim

	options_group NLPTester {
        print_all = false;
        throw_exception = true;
	}
 \endverbatim
 */
class NLPTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPTester >
{
public:

	///
	NLPTesterSetOptions(
		  NLPTester* target = 0
		, const char opt_grp_name[] = "NLPTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class NLPTesterSetOptions

}	// end namespace NLPInterfacePack

#endif	// NLP_TESTER_SET_OPTIONS_H
