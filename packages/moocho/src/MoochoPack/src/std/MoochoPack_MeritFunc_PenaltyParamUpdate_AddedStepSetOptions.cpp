// ////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.cpp
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

#include "../std/MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 4;

	const char options_group_name[] = "MeritFuncPenaltyParamUpdate";

	enum local_EOptions {
	    SMALL_MU,
	    MIN_MU_RATIO,
	    MULT_FACTOR,
	    KKT_NEAR_SOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "small_mu",
	    "min_mu_ratio",
	    "mult_factor",
	    "kkt_near_sol"
	};

}

namespace MoochoPack {

MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
			MeritFunc_PenaltyParamUpdate_AddedStep* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			MeritFunc_PenaltyParamUpdate_AddedStep >( target )
{}

void MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case SMALL_MU: {
			target().small_mu( ::atof( option_value.c_str() ) );
			break;
		}
		case MIN_MU_RATIO: {
			target().min_mu_ratio( ::atof( option_value.c_str() ) );
			break;
		}
		case MULT_FACTOR: {
			target().mult_factor( ::atof( option_value.c_str() ) );
			break;
		}
		case KKT_NEAR_SOL: {
			target().kkt_near_sol( ::atof( option_value.c_str() ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
