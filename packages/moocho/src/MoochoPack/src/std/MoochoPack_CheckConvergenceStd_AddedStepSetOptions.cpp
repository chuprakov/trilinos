// ////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include "../../include/std/CheckConvergenceStd_AddedStepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 1;

	enum local_EOptions {
		SCALE_KKT_ERROR_BY
	};

	const char* local_SOptions[local_num_options]	= {
		"scale_kkt_error_by",
	};

}

namespace ReducedSpaceSQPPack {

CheckConvergenceStd_AddedStepSetOptions::CheckConvergenceStd_AddedStepSetOptions(
			  CheckConvergenceStd_AddedStep* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			CheckConvergenceStd_AddedStep >( target )
{}

void CheckConvergenceStd_AddedStepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef CheckConvergenceStd_AddedStep target_t;
	switch( (local_EOptions)option_num ) {
	    case SCALE_KKT_ERROR_BY:
		{
			const std::string &option = option_value.c_str();
			if( option == "SCALE_BY_ONE" )
				target().scale_kkt_error_by( target_t::SCALE_BY_ONE );
			else if( option == "SCALE_BY_NORM_2_X" )
				target().scale_kkt_error_by( target_t::SCALE_BY_NORM_2_X );
			else if( option == "SCALE_BY_NORM_INF_X" )
				target().scale_kkt_error_by( target_t::SCALE_BY_NORM_INF_X );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"scale_kkt_error_by\".  Only the options "
					"SCALE_BY_ONE, SCALE_BY_NORM_2_X, and SCALE_BY_NORM_INF_X "
					"are available" );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 