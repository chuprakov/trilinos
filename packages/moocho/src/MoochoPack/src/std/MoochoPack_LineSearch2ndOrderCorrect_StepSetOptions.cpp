// ////////////////////////////////////////////////////////////////
// LineSearch2ndOrderCorrect_StepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <math.h>

#include "../../include/std/LineSearch2ndOrderCorrect_StepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 8;

	enum local_EOptions {
	    NEWTON_OLEVEL,
		CONSTR_NORM_THRESHOLD,
		CONSTR_INCR_RATIO,
		AFTER_K_ITER,
		FORCED_CONSTR_REDUCTION,
		FORCED_REDUCT_RATIO,
		MAX_STEP_RATIO,
		MAX_NEWTON_ITER
	};

	const char* local_SOptions[local_num_options]	= {
	    "newton_olevel",
		"constr_norm_threshold",
		"constr_incr_ratio",
		"after_k_iter",
		"forced_constr_reduction",
		"forced_reduct_ratio",
		"max_step_ratio",
		"max_newton_iter"
	};

}

namespace ReducedSpaceSQPPack {

LineSearch2ndOrderCorrect_StepSetOptions::LineSearch2ndOrderCorrect_StepSetOptions(
			  LineSearch2ndOrderCorrect_Step* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			LineSearch2ndOrderCorrect_Step >( target )
{}

void LineSearch2ndOrderCorrect_StepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	typedef LineSearch2ndOrderCorrect_Step target_t;
	switch( (local_EOptions)option_num ) {
	    case NEWTON_OLEVEL:
		{
			const std::string &option = option_value.c_str();
			if( option == "PRINT_USE_DEFAULT" )
				target().newton_olevel( target_t::PRINT_USE_DEFAULT );
			else if( option == "PRINT_NOTHING" )
				target().newton_olevel( target_t::PRINT_NEWTON_NOTHING );
			else if( option == "PRINT_SUMMARY_INFO" )
				target().newton_olevel( target_t::PRINT_NEWTON_SUMMARY_INFO );
			else if( option == "PRINT_STEPS" )
				target().newton_olevel( target_t::PRINT_NEWTON_STEPS );
			else if( option == "PRINT_VECTORS" )
				target().newton_olevel( target_t::PRINT_NEWTON_VECTORS );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"newton_olevel\"." );
			break;
		}
	    case CONSTR_NORM_THRESHOLD:
			target().constr_norm_threshold(::fabs(::atof(option_value.c_str())));
			break;
	    case CONSTR_INCR_RATIO:
			target().constr_incr_ratio(::fabs(::atof(option_value.c_str())));
			break;
		case AFTER_K_ITER:
			target().after_k_iter(::abs(::atoi(option_value.c_str())));
			break;
		case FORCED_CONSTR_REDUCTION:
		{
			const std::string &option = option_value.c_str();
			if( option == "LESS_X_D" )
				target().forced_constr_reduction(target_t::CONSTR_LESS_X_D );
			else if( option == "LESS_X" )
				target().forced_constr_reduction( target_t::CONSTR_LESS_X );
			else
				throw std::invalid_argument( "Error, incorrect value for "
					"\"forced_constr_reduction\"." );
			break;
		}
	    case FORCED_REDUCT_RATIO:
			target().forced_reduct_ratio(::fabs(::atof(option_value.c_str())));
			break;
		case MAX_STEP_RATIO:
			target().max_step_ratio(::fabs(::atof(option_value.c_str())));
			break;
		case MAX_NEWTON_ITER:
			target().max_newton_iter(::abs(::atoi(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 
