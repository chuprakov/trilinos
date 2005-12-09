// ////////////////////////////////////////////////////////////////
// LineSearchWatchDog_StepSetOptions.cpp
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

#include "../std/MoochoPack_LineSearchWatchDog_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 2;

	const char options_group_name[] = "LineSearchWatchDog";

	enum local_EOptions {
		OPT_KKT_ERR_THESHOLD,
		FEAS_KKT_ERR_THESHOLD
	};

	const char* local_SOptions[local_num_options]	= {
		"opt_kkt_err_threshold",
		"feas_kkt_err_threshold"
	};

}

namespace MoochoPack {

LineSearchWatchDog_StepSetOptions::LineSearchWatchDog_StepSetOptions(
			LineSearchWatchDog_Step* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			LineSearchWatchDog_Step >( target )
{}

void LineSearchWatchDog_StepSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case OPT_KKT_ERR_THESHOLD: {
			target().opt_kkt_err_threshold( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case FEAS_KKT_ERR_THESHOLD: {
			target().feas_kkt_err_threshold( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
