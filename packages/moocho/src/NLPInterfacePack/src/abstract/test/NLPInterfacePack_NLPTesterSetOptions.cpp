// ////////////////////////////////////////////////////////////////
// NLPTesterSetOptions.cpp
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

#include "NLPInterfacePack_NLPTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 2;

	enum local_EOptions {
		PRINT_ALL
		,TEST_FOR_EXCEPTION
	};

	const char* local_SOptions[local_num_options]	= {
		"print_all"
		,"throw_exception"
	};

}

namespace NLPInterfacePack {

NLPTesterSetOptions::NLPTesterSetOptions(
			  NLPTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPTester >( target )
{}

void NLPTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::StringToBool;

	typedef NLPTester target_t;
	switch( (local_EOptions)option_num ) {
		case PRINT_ALL:
			target().print_all(
				StringToBool( "print_all", option_value.c_str() )
				);
			break;
		case TEST_FOR_EXCEPTION:
			target().throw_exception(
				StringToBool( "throw_exception", option_value.c_str() )
				);
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace NLPInterfacePack
