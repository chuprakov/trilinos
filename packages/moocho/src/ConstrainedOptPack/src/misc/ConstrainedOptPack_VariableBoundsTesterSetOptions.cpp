// ////////////////////////////////////////////////////////////////
// VariableBoundsTesterSetOptions.cpp
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

#include "ConstrainedOptPack_VariableBoundsTesterSetOptions.hpp"

// Define the options
namespace {

	const int local_num_options = 2;

	enum local_EOptions {
		WARNING_TOL
		,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "warning_tol"
	    ,"error_tol"
	};

}

namespace ConstrainedOptPack {

VariableBoundsTesterSetOptions::VariableBoundsTesterSetOptions(
			  VariableBoundsTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			VariableBoundsTester >( target )
{}

void VariableBoundsTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	typedef VariableBoundsTester target_t;
	switch( (local_EOptions)option_num ) {
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

}	// end namespace ConstrainedOptPack
