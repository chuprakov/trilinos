// ////////////////////////////////////////////////////////////////
// LineSearchWatchDog_StepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "../../include/std/LineSearchWatchDog_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

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

namespace ReducedSpaceSQPPack {

LineSearchWatchDog_StepSetOptions::LineSearchWatchDog_StepSetOptions(
			LineSearchWatchDog_Step* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			LineSearchWatchDog_Step >( target )
{}

void LineSearchWatchDog_StepSetOptions::set_option(
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

}	// end namespace ReducedSpaceSQPPack 
