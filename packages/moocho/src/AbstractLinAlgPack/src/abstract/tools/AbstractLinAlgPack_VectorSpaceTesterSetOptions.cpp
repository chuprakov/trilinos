// ////////////////////////////////////////////////////////////////
// VectorSpaceTesterSetOptions.cpp
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

#include <assert.h>
#include <math.h>

#include "AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 6;

	enum local_EOptions {
		PRINT_ALL_TESTS
		,PRINT_VECTORS
		,TEST_FOR_EXCEPTION
		,NUM_RANDOM_TESTS
	    ,WARNING_TOL
	    ,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"print_all_tests"
		,"print_vectors"
		,"throw_exception"
		,"num_random_tests"
	    ,"warning_tol"
	    ,"error_tol"
	};

}

namespace AbstractLinAlgPack {

VectorSpaceTesterSetOptions::VectorSpaceTesterSetOptions(
			  VectorSpaceTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			VectorSpaceTester >( target )
{}

void VectorSpaceTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;
	typedef VectorSpaceTester target_t;
	switch( (local_EOptions)option_num ) {
		case PRINT_ALL_TESTS:
			target().print_all_tests(
				StringToBool( "print_all_tests", option_value.c_str() )
				);
			break;
		case PRINT_VECTORS:
			target().print_vectors(
				StringToBool( "print_vectors", option_value.c_str() )
				);
			break;
		case TEST_FOR_EXCEPTION:
			target().throw_exception(
				StringToBool( "throw_exception", option_value.c_str() )
				);
			break;
	    case NUM_RANDOM_TESTS:
			target().num_random_tests(::abs(::atoi(option_value.c_str())));
			break;
	    case WARNING_TOL:
			target().warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case ERROR_TOL:
			target().error_tol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace AbstractLinAlgPack
