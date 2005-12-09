// ////////////////////////////////////////////////////////////////
// QPSolverRelaxedTesterSetOptions.cpp
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

#include "ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.hpp"

// Define the options
namespace {

	const int local_num_options = 6;

	enum local_EOptions {
		OPT_WARNING_TOL
		,OPT_ERROR_TOL
		,FEAS_WARNING_TOL
		,FEAS_ERROR_TOL
		,COMP_WARNING_TOL
		,COMP_ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "opt_warning_tol"
	    ,"opt_error_tol"
	    ,"feas_warning_tol"
	    ,"feas_error_tol"
	    ,"comp_warning_tol"
	    ,"comp_error_tol"
	};

}

namespace ConstrainedOptPack {

QPSolverRelaxedTesterSetOptions::QPSolverRelaxedTesterSetOptions(
			  QPSolverRelaxedTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			QPSolverRelaxedTester >( target )
{}

void QPSolverRelaxedTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	typedef QPSolverRelaxedTester target_t;
	switch( (local_EOptions)option_num ) {
	    case OPT_WARNING_TOL:
			target().opt_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case OPT_ERROR_TOL:
			target().opt_error_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case FEAS_WARNING_TOL:
			target().feas_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case FEAS_ERROR_TOL:
			target().feas_error_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case COMP_WARNING_TOL:
			target().comp_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case COMP_ERROR_TOL:
			target().comp_error_tol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ConstrainedOptPack
