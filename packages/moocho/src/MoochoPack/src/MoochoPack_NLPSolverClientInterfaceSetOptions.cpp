// ////////////////////////////////////////////////////////////////
// rSQPSolverClientInterfaceSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include "../include/rSQPSolverClientInterfaceSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 6;

	enum local_EOptions {
        MAX_ITER,
        MAX_RUN_TIME,
        OPT_TOL,
        FEAS_TOL,
        STEP_TOL,
        MAX_VAR_BOUNDS_VIOL
	};

	const char* local_SOptions[local_num_options]	= {
        "max_iter",
        "max_run_time",
        "opt_tol",
        "feas_tol",
        "step_tol",
        "max_var_bounds_viol"
	};

}

namespace ReducedSpaceSQPPack {

rSQPSolverClientInterfaceSetOptions::rSQPSolverClientInterfaceSetOptions(
			  rSQPSolverClientInterface* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			rSQPSolverClientInterface >( target )
{}

void rSQPSolverClientInterfaceSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef rSQPSolverClientInterface target_t;
	switch( (local_EOptions)option_num ) {
	    case MAX_ITER:
			target().max_iter(::abs(::atoi(option_value.c_str())));
			break;
		case MAX_RUN_TIME:
			target().max_run_time(::fabs(::atof(option_value.c_str())));
			break;
		case OPT_TOL:
			target().opt_tol(::fabs(::atof(option_value.c_str())));
			break;
		case FEAS_TOL:
			target().feas_tol(::fabs(::atof(option_value.c_str())));
			break;
		case STEP_TOL:
			target().step_tol(::fabs(::atof(option_value.c_str())));
			break;
		case MAX_VAR_BOUNDS_VIOL:
			target().max_var_bounds_viol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 