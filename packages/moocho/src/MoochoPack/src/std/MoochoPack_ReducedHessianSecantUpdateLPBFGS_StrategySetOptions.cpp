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

	const int local_num_options = 7;

	enum local_EOptions {
		ACT_SET_FRAC_PROJ_START
        ,SUPER_BASIC_MULT_DROP_TOL
		,NUM_SUPERBASICS_SWITCH_DENSE
		,ACT_SET_FRAC_SWITCH_DENSE
		,MIN_NUM_UPDATES_SWITCH_DENSE
		,MAX_NUM_UPDATES_SWITCH_DENSE
		,NUM_ADD_RECENT_UPDATES
	};

	const char* local_SOptions[local_num_options]	= {
		"act_set_frac_proj_start"
        ,"super_basic_mult_drop_tol"
		,"num_superbasics_switch_dense"
		,"act_set_frac_switch_dense"
		,"min_num_updates_switch_dense"
		,"max_num_updates_switch_dense"
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
		case ACT_SET_FRAC_PROJ_START: {
			target().act_set_frac_proj_start( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case SUPER_BASIC_MULT_DROP_TOL: {
			target().super_basic_mult_drop_tol( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case NUM_SUPERBASICS_SWITCH_DENSE: {
			target().num_superbasics_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		case ACT_SET_FRAC_SWITCH_DENSE: {
			target().act_set_frac_switch_dense( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case MIN_NUM_UPDATES_SWITCH_DENSE: {
			target().min_num_updates_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
			break;
		}
		case MAX_NUM_UPDATES_SWITCH_DENSE: {
			target().max_num_updates_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
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
