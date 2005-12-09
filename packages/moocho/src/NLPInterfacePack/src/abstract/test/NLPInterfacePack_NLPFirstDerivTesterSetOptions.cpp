// ////////////////////////////////////////////////////////////////
// NLPFirstDerivTesterSetOptions.cpp
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

#include "NLPInterfacePack_NLPFirstDerivTesterSetOptions.hpp"

// Define the options
namespace {

	const int local_num_options = 4;

	enum local_EOptions {
		FD_TESTING_METHOD
		,NUM_FD_DIRECTIONS
	    ,WARNING_TOL
	    ,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "fd_testing_method"
	    ,"num_fd_directions"
	    ,"warning_tol"
	    ,"error_tol"
	};

}

namespace NLPInterfacePack {

NLPFirstDerivTesterSetOptions::NLPFirstDerivTesterSetOptions(
			  NLPFirstDerivTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPFirstDerivTester >( target )
{}

void NLPFirstDerivTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	typedef NLPFirstDerivTester target_t;
	switch( (local_EOptions)option_num ) {
	    case FD_TESTING_METHOD:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_COMPUTE_ALL" )
				target().fd_testing_method( target_t::FD_COMPUTE_ALL );
			else if( option == "FD_DIRECTIONAL" )
				target().fd_testing_method( target_t::FD_DIRECTIONAL );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"fd_testing_method\".  Only the options "
					"FD_COMPUTE_ALL and FD_DIRECTIONAL are available" );
			break;
		}
	    case NUM_FD_DIRECTIONS:
			target().num_fd_directions(::abs(::atoi(option_value.c_str())));
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

}	// end namespace NLPInterfacePack
