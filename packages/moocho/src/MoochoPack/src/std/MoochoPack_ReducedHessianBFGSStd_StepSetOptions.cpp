// ////////////////////////////////////////////////////////////////
// ReducedHessianBFGSStd_StepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "../../include/std/ReducedHessianBFGSStd_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 5;

	enum local_EOptions {
		RESCALE_INIT_IDENTITY
		,USE_DAMPENING
		,SECANT_TESTING
		,SECANT_WARNING_TOL
		,SECANT_ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"rescale_init_identity"
	    ,"use_dampening"
		,"secant_testing"
		,"secant_warning_tol"
	    ,"secant_error_tol"
	};

}

namespace ReducedSpaceSQPPack {

ReducedHessianBFGSStd_StepSetOptions::ReducedHessianBFGSStd_StepSetOptions(
			  ReducedHessianBFGSStd_Step* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			ReducedHessianBFGSStd_Step >( target )
{}

void ReducedHessianBFGSStd_StepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;
	typedef ReducedHessianBFGSStd_Step target_t;
	switch( (local_EOptions)option_num ) {
	    case RESCALE_INIT_IDENTITY:
			target().rescale_init_identity(
				StringToBool( "rescale_init_identity", option_value.c_str() ));
			break;
	    case USE_DAMPENING:
			target().use_dampening(
				StringToBool( "use_dampening", option_value.c_str() ));
			break;
	    case SECANT_TESTING:
		{
			const std::string &option = option_value.c_str();
			if( option == "DEFAULT" )
				target().secant_testing( target_t::SECANT_TEST_DEFAULT );
			else if( option == "TEST" )
				target().secant_testing( target_t::SECANT_TEST_ALWAYS );
			else if( option == "NO_TEST" )
				target().secant_testing( target_t::SECANT_NO_TEST );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"secant_testing\"." );
			break;
		}
	    case SECANT_WARNING_TOL:
			target().secant_warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case SECANT_ERROR_TOL:
			target().secant_error_tol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 