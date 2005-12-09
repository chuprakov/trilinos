// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.cpp
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

#include "../std/MoochoPack_ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 3;

	enum local_EOptions {
		ACT_SET_FRAC_PROJ_START
		,PROJECT_ERROR_TOL
		,SUPER_BASIC_MULT_DROP_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"act_set_frac_proj_start"
		,"project_error_tol"
		,"super_basic_mult_drop_tol"
	};

}

namespace MoochoPack {

ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions::ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions(
	ReducedHessianSecantUpdateBFGSProjected_Strategy* target
	, const char opt_grp_name[] )
	: OptionsFromStreamPack::SetOptionsFromStreamNode(
		opt_grp_name, local_num_options, local_SOptions )
	, OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateBFGSProjected_Strategy >( target )
{}

void ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions::setOption(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case ACT_SET_FRAC_PROJ_START: {
			target().act_set_frac_proj_start( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case PROJECT_ERROR_TOL: {
			target().project_error_tol( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case SUPER_BASIC_MULT_DROP_TOL: {
			target().super_basic_mult_drop_tol( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
