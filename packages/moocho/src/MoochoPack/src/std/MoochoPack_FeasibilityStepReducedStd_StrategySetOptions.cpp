// ////////////////////////////////////////////////////////////////
// FeasibilityStepReducedStd_StrategySetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>
#include <math.h>

#include "../../include/std/FeasibilityStepReducedStd_StrategySetOptions.h"

// Define the options
namespace {

	const int local_num_options = 2;

	enum local_EOptions {
		QP_OBJECTIVE
		,QP_TESTING
	};

	const char* local_SOptions[local_num_options]	= {
		"qp_objective"
		,"qp_testing"
	};

}

namespace ReducedSpaceSQPPack {

FeasibilityStepReducedStd_StrategySetOptions::FeasibilityStepReducedStd_StrategySetOptions(
			  FeasibilityStepReducedStd_Strategy* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			FeasibilityStepReducedStd_Strategy >( target )
{}

void FeasibilityStepReducedStd_StrategySetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef FeasibilityStepReducedStd_Strategy target_t;
	switch( (local_EOptions)option_num ) {
	    case QP_OBJECTIVE:
		{
			const std::string &option = option_value.c_str();
			if( option == "OBJ_MIN_FULL_STEP" )
				target().qp_objective( target_t::OBJ_MIN_FULL_STEP );
			else if( option == "OBJ_MIN_NULL_SPACE_STEP" )
				target().qp_objective( target_t::OBJ_MIN_NULL_SPACE_STEP );
			else if( option == "OBJ_RSQP" )
				target().qp_objective( target_t::OBJ_RSQP );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"qp_objective\".  Only the options "
					"OBJ_MIN_FULL_STEP, OBJ_MIN_NULL_SPACE_STEP, and  OBJ_RSQP"
					"are available" );
			break;
		}
	    case QP_TESTING:
		{
			const std::string &option = option_value.c_str();
			if( option == "QP_TEST_DEFAULT" )
				target().qp_testing( target_t::QP_TEST_DEFAULT );
			else if( option == "QP_TEST" )
				target().qp_testing( target_t::QP_TEST );
			else if( option == "QP_NO_TEST" )
				target().qp_testing( target_t::QP_NO_TEST );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"qp_testing\".  Only the options "
					"QP_TEST_DEFAULT, QP_TEST, and QP_NO_TEST "
					"are available" );
			break;
		}
	    default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 
