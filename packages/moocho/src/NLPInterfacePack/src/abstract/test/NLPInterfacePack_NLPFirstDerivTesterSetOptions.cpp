// ////////////////////////////////////////////////////////////////
// NLPFirstDerivativesTesterSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "NLPFirstDerivativesTesterSetOptions.h"

// Define the options
namespace {

	const int local_num_options = 4;

	enum local_EOptions {
		FD_TESTING_METHOD
		,NUM_FD_DIRECTIONS
	    ,WARNING_TOL
	    ,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "fd_testing_method"
	    ,"num_fd_directions"
	    ,"warning_tol"
	    ,"error_tol"
	};

}

namespace NLPInterfacePack {
namespace TestingPack {

NLPFirstDerivativesTesterSetOptions::NLPFirstDerivativesTesterSetOptions(
			  NLPFirstDerivativesTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPFirstDerivativesTester >( target )
{}

void NLPFirstDerivativesTesterSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef NLPFirstDerivativesTester target_t;
	switch( (local_EOptions)option_num ) {
	    case FD_TESTING_METHOD:
		{
			const std::string &option = option_value.c_str();
			if( option == "COMPUTE_ALL" )
				target().fd_testing_method( target_t::FD_COMPUTE_ALL );
			else if( option == "DIRECTIONAL" )
				target().fd_testing_method( target_t::FD_DIRECTIONAL );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"fd_testing_method\".  Only the options "
					"COMPUTE_ALL and DIRECTIONAL are available" );
			break;
		}
	    case NUM_FD_DIRECTIONS:
			target().num_fd_directions(::fabs(::atoi(option_value.c_str())));
			break;
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

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack
