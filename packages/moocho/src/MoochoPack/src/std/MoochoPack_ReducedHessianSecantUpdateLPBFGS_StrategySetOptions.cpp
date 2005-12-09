// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.cpp
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

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 4;

	enum local_EOptions {
		MIN_NUM_UPDATES_PROJ_START
		,MAX_NUM_UPDATES_PROJ_START
		,NUM_SUPERBASICS_SWITCH_DENSE
		,NUM_ADD_RECENT_UPDATES
	};

	const char* local_SOptions[local_num_options]	= {
		"min_num_updates_proj_start"
		,"max_num_updates_proj_start"
		,"num_superbasics_switch_dense"
		,"num_add_recent_updates"
	};

}

namespace MoochoPack {

ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
	ReducedHessianSecantUpdateLPBFGS_Strategy* target
	, const char opt_grp_name[] )
	: OptionsFromStreamPack::SetOptionsFromStreamNode(
		opt_grp_name, local_num_options, local_SOptions )
	, OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateLPBFGS_Strategy >( target )
{}

void ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::setOption(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case MIN_NUM_UPDATES_PROJ_START: {
			target().min_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		case MAX_NUM_UPDATES_PROJ_START: {
			target().max_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		case NUM_SUPERBASICS_SWITCH_DENSE: {
			target().num_superbasics_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		case NUM_ADD_RECENT_UPDATES: {
			target().num_add_recent_updates( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
