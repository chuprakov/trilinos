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

	const int local_num_options = 3;

	enum local_EOptions {
		ACT_SET_FRAC_PROJ_START
		,PROJECT_ERROR_TOL
		,SUPER_BASIC_MULT_DROP_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"act_set_frac_proj_start"
		,"project_error_tol"
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
		case ACT_SET_FRAC_PROJ_START: {
			target().act_set_frac_proj_start( ::fabs( ::atof( option_value.c_str() ) ) );
			break;
		}
		case PROJECT_ERROR_TOL: {
			target().project_error_tol( ::fabs( ::atof( option_value.c_str() ) ) );
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
