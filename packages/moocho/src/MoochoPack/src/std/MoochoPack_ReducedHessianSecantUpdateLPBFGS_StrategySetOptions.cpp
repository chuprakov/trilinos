// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.h"
#include "Misc/include/StringToBool.h"

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

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
	ReducedHessianSecantUpdateLPBFGS_Strategy* target
	, const char opt_grp_name[] )
	: OptionsFromStreamPack::SetOptionsFromStreamNode(
		opt_grp_name, local_num_options, local_SOptions )
	, OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateLPBFGS_Strategy >( target )
{}

void ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::set_option(
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

}	// end namespace ReducedSpaceSQPPack 
