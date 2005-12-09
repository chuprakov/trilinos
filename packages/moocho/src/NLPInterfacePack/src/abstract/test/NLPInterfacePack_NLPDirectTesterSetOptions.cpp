// ////////////////////////////////////////////////////////////////
// NLPDirectTesterSetOptions.cpp
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

#include "NLPInterfacePack_NLPDirectTesterSetOptions.hpp"

// Define the options
namespace {

	const int local_num_options = 7;

	enum local_EOptions {
		GF_TESTING_METHOD
	    ,GF_WARNING_TOL
	    ,GF_ERROR_TOL
		,GC_TESTING_METHOD
	    ,GC_WARNING_TOL
	    ,GC_ERROR_TOL
		,NUM_FD_DIRECTIONS
	};

	const char* local_SOptions[local_num_options]	= {
	    "Gf_testing_method"
	    ,"Gf_warning_tol"
	    ,"Gf_error_tol"
	    ,"Gc_testing_method"
	    ,"Gc_warning_tol"
	    ,"Gc_error_tol"
	    ,"num_fd_directions"
	};

}

namespace NLPInterfacePack {

NLPDirectTesterSetOptions::NLPDirectTesterSetOptions(
			  NLPDirectTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPDirectTester >( target )
{}

void NLPDirectTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	typedef NLPDirectTester target_t;
	switch( (local_EOptions)option_num ) {
	    case GF_TESTING_METHOD:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_COMPUTE_ALL" )
				target().Gf_testing_method( target_t::FD_COMPUTE_ALL );
			else if( option == "FD_DIRECTIONAL" )
				target().Gf_testing_method( target_t::FD_DIRECTIONAL );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"Gf_testing_method\".  Only the options "
					"FD_COMPUTE_ALL and FD_DIRECTIONAL are available" );
			break;
		}
	    case GF_WARNING_TOL:
			target().Gf_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case GF_ERROR_TOL:
			target().Gf_error_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case GC_TESTING_METHOD:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_COMPUTE_ALL" )
				target().Gc_testing_method( target_t::FD_COMPUTE_ALL );
			else if( option == "FD_DIRECTIONAL" )
				target().Gc_testing_method( target_t::FD_DIRECTIONAL );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"Gc_testing_method\".  Only the options "
					"FD_COMPUTE_ALL and FD_DIRECTIONAL are available" );
			break;
		}
	    case GC_WARNING_TOL:
			target().Gc_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case GC_ERROR_TOL:
			target().Gc_error_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case NUM_FD_DIRECTIONS:
			target().num_fd_directions(::abs(::atoi(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace NLPInterfacePack
