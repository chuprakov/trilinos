// ////////////////////////////////////////////////////////////////
// EvalNewPointStd_StepSetOptions.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include "../../include/std/EvalNewPointStd_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 1;

	enum local_EOptions {
		FD_DERIV_TESTING
	};

	const char* local_SOptions[local_num_options]	= {
		"fd_deriv_testing"
	};

}

namespace ReducedSpaceSQPPack {

EvalNewPointStd_StepSetOptions::EvalNewPointStd_StepSetOptions(
			  EvalNewPointStd_Step* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			EvalNewPointStd_Step >( target )
{}

void EvalNewPointStd_StepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;

	typedef EvalNewPointStd_Step target_t;
	switch( (local_EOptions)option_num ) {
	    case FD_DERIV_TESTING:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_DEFAULT" )
				target().fd_deriv_testing( target_t::FD_DEFAULT );
			else if( option == "FD_TEST" )
				target().fd_deriv_testing( target_t::FD_TEST );
			else if( option == "FD_NO_TEST" )
				target().fd_deriv_testing( target_t::FD_NO_TEST );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"fd_deriv_testing\".  Only the options "
					"FD_DEFAULT, FD_TEST, and FD_NO_TEST "
					"are available" );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 
