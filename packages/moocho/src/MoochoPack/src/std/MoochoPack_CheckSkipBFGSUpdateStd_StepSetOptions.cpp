// ////////////////////////////////////////////////////////////////
// CheckSkipBFGSUpdateStd_StepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "../../include/std/CheckSkipBFGSUpdateStd_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

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

namespace ReducedSpaceSQPPack {

CheckSkipBFGSUpdateStd_StepSetOptions::CheckSkipBFGSUpdateStd_StepSetOptions(
			CheckSkipBFGSUpdateStd_Step* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			CheckSkipBFGSUpdateStd_Step >( target )
{}

void CheckSkipBFGSUpdateStd_StepSetOptions::set_option(
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

}	// end namespace ReducedSpaceSQPPack 