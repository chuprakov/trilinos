// ////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp
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

#ifndef VECTOR_SPACE_TESTER_SET_OPTIONS_H
#define VECTOR_SPACE_TESTER_SET_OPTIONS_H

#include "AbstractLinAlgPack_VectorSpaceTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace AbstractLinAlgPack {

///
/** Set options for VectorSpaceTester from an
  * OptionsFromStream object.
  *
  * The default options group name is VectorSpaceTester.
  *
  * The options group is:
  *
  \verbatim

    options_group VectorSpaceTester {
    *    print_all_tests = true;
        print_all_tests = false;
    *    print_vectors = true;
        print_vectors = false;
        throw_exception = true;
    *    throw_exception = false;
        num_random_tests = 4;
        warning_tol = 1e-14;
        error_tol   = 1e-10;
    }
  \endverbatim
  */
class VectorSpaceTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			VectorSpaceTester >
{
public:

	///
	VectorSpaceTesterSetOptions(
		  VectorSpaceTester* target = 0
		, const char opt_grp_name[] = "VectorSpaceTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class VectorSpaceTesterSetOptions

}	// end namespace AbstractLinAlgPack

#endif	// VECTOR_SPACE_TESTER_SET_OPTIONS_H
