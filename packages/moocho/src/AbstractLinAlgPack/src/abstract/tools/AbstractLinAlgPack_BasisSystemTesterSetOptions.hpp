// ////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp
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

#ifndef BASIS_SYSTEM_TESTER_SET_OPTIONS_H
#define BASIS_SYSTEM_TESTER_SET_OPTIONS_H

#include "AbstractLinAlgPack_BasisSystemTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace AbstractLinAlgPack {

///
/** Set options for BasisSystemTester from an
  * OptionsFromStream object.
  *
  * The default options group name is BasisSystemTester.
  *
  * The options group is:
  *
  \verbatim

    options_group BasisSystemTester {
        print_tests = PRINT_NONE;
    *    print_tests = PRINT_BASIC;
    *    print_tests = PRINT_MORE;
    *    print_tests = PRINT_ALL;
    *    dump_all = true;
        dump_all = false;
        throw_exception = true;
    *    throw_exception = false;
        num_random_tests = 1;
        warning_tol = 1e-14;
        error_tol   = 1e-10;
    }
  \endverbatim
  */
class BasisSystemTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			BasisSystemTester >
{
public:

	///
	BasisSystemTesterSetOptions(
		  BasisSystemTester* target = 0
		, const char opt_grp_name[] = "BasisSystemTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class BasisSystemTesterSetOptions

}	// end namespace AbstractLinAlgPack

#endif	// BASIS_SYSTEM_TESTER_SET_OPTIONS_H
