// ////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_StepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 1;

	enum local_EOptions {
		FD_DERIV_TESTING
	};

	const char* local_SOptions[local_num_options]	= {
		"fd_deriv_testing"
	};

}

namespace ReducedSpaceSQPPack {

EvalNewPointTailoredApproach_StepSetOptions::EvalNewPointTailoredApproach_StepSetOptions(
			  EvalNewPointTailoredApproach_Step* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			EvalNewPointTailoredApproach_Step >( target )
{}

void EvalNewPointTailoredApproach_StepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;

	typedef EvalNewPointTailoredApproach_Step target_t;
	switch( (local_EOptions)option_num ) {
	    case FD_DERIV_TESTING:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_DEFAULT" )
				target().fd_deriv_testing( target_t::FD_DEFAULT );
			else if( option == "FD_TEST" )
				target().fd_deriv_testing( target_t::FD_TEST );
			else if( option == "FD_NO_TEST" )
				target().fd_deriv_testing( target_t::FD_NO_TEST );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"fd_deriv_testing\".  Only the options "
					"FD_DEFAULT, FD_TEST, and FD_NO_TEST "
					"are available" );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 