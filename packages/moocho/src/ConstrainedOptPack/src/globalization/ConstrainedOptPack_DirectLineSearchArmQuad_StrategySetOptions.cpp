// ////////////////////////////////////////////////////////////////
// DirectLineSearchArmQuad_StrategySetOptions.cpp
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

#include "DirectLineSearchArmQuad_StrategySetOptions.h"
#include "StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 4;

	enum local_EOptions {
		SLOPE_FRAC,
		MIN_FRAC_STEP,
		MAX_FRAC_STEP,
		MAX_LS_ITER
	};

	const char* local_SOptions[local_num_options]	= {
		"slope_frac",
		"min_frac_step",
		"max_frac_step",
		"max_ls_iter"
	};

}

namespace ConstrainedOptimizationPack {

DirectLineSearchArmQuad_StrategySetOptions::DirectLineSearchArmQuad_StrategySetOptions(
			  DirectLineSearchArmQuad_Strategy* qp_solver
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			DirectLineSearchArmQuad_Strategy >( qp_solver )
{}

void DirectLineSearchArmQuad_StrategySetOptions::set_option(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case SLOPE_FRAC: {
			target().eta( ::atof( option_value.c_str() ) );
			break;
		}
		case MIN_FRAC_STEP: {
			target().min_frac( ::atof( option_value.c_str() ) );
			break;
		}
		case MAX_FRAC_STEP: {
			target().max_frac( ::atof( option_value.c_str() ) );
			break;
		}
		case MAX_LS_ITER: {
			target().set_max_iter( ::atof( option_value.c_str() ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ConstrainedOptimizationPack 
