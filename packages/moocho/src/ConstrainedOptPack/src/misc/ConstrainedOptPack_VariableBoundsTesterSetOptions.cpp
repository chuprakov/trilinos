// ////////////////////////////////////////////////////////////////
// VariableBoundsTesterSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "ConstrainedOptimizationPack/include/VariableBoundsTesterSetOptions.h"

// Define the options
namespace {

	const int local_num_options = 2;

	enum local_EOptions {
		WARNING_TOL
		,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "warning_tol"
	    ,"error_tol"
	};

}

namespace ConstrainedOptimizationPack {

VariableBoundsTesterSetOptions::VariableBoundsTesterSetOptions(
			  VariableBoundsTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			VariableBoundsTester >( target )
{}

void VariableBoundsTesterSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef VariableBoundsTester target_t;
	switch( (local_EOptions)option_num ) {
	    case WARNING_TOL:
			target().warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case ERROR_TOL:
			target().error_tol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ConstrainedOptimizationPack
