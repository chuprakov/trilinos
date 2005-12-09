// ////////////////////////////////////////////////////////////////
// CheckSkipBFGSUpdateStd_StepSetOptions.cpp
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

#include "../std/MoochoPack_CheckSkipBFGSUpdateStd_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 1;

	const char options_group_name[] = "CheckSkipBFGSUpdateStd";

	enum local_EOptions {
		SKIP_BFGS_PROP_CONST
	};

	const char* local_SOptions[local_num_options]	= {
		"skip_bfgs_prop_const"
	};

}

namespace MoochoPack {

CheckSkipBFGSUpdateStd_StepSetOptions::CheckSkipBFGSUpdateStd_StepSetOptions(
			CheckSkipBFGSUpdateStd_Step* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			CheckSkipBFGSUpdateStd_Step >( target )
{}

void CheckSkipBFGSUpdateStd_StepSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case SKIP_BFGS_PROP_CONST: {
			target().skip_bfgs_prop_const( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
