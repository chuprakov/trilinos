// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "../../include/std/ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 2;

	enum local_EOptions {
		PROJ_START_ACT_SET_FRAC
		,SUPER_BASIC_MULT_DROP_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"proj_start_act_set_frac"
		,"super_basic_mult_drop_tol"
	};

}

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions::ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions(
	ReducedHessianSecantUpdateBFGSProjected_Strategy* target
	, const char opt_grp_name[] )
	: OptionsFromStreamPack::SetOptionsFromStreamNode(
		opt_grp_name, local_num_options, local_SOptions )
	, OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateBFGSProjected_Strategy >( target )
{}

void ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions::set_option(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case PROJ_START_ACT_SET_FRAC: {
			target().proj_start_act_set_frac( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case SUPER_BASIC_MULT_DROP_TOL: {
			target().super_basic_mult_drop_tol( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 
