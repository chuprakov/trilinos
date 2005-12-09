// ////////////////////////////////////////////////////////////////
// CalcFiniteDiffProdSetOptions.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <assert.h>
#include <math.h>

#include "NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

	const int local_num_options = 6;

	enum local_EOptions {
		FD_METHOD_ORDER
		,FD_STEP_SELECT
		,FD_STEP_SIZE
		,FD_STEP_SIZE_MIN
		,FD_STEP_SIZE_F
		,FD_STEP_SIZE_C
	};

	const char* local_SOptions[local_num_options]	= {
		"fd_method_order"
		,"fd_step_select"
		,"fd_step_size"
		,"fd_step_size_min"
		,"fd_step_size_f"
		,"fd_step_size_c"
	};

}

namespace NLPInterfacePack {

CalcFiniteDiffProdSetOptions::CalcFiniteDiffProdSetOptions(
	CalcFiniteDiffProd* target
	,const char opt_grp_name[]
	)
	:OptionsFromStreamPack::SetOptionsFromStreamNode(opt_grp_name,local_num_options,local_SOptions)
	,OptionsFromStreamPack::SetOptionsToTargetBase<CalcFiniteDiffProd>(target)
{}

void CalcFiniteDiffProdSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	typedef CalcFiniteDiffProd target_t;
	switch( (local_EOptions)option_num ) {
	    case FD_METHOD_ORDER:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_ORDER_ONE" )
				target().fd_method_order( target_t::FD_ORDER_ONE );
			else if( option == "FD_ORDER_TWO" )
				target().fd_method_order( target_t::FD_ORDER_TWO );
			else if( option == "FD_ORDER_TWO_CENTRAL" )
				target().fd_method_order( target_t::FD_ORDER_TWO_CENTRAL );
			else if( option == "FD_ORDER_TWO_AUTO" )
				target().fd_method_order( target_t::FD_ORDER_TWO_AUTO );
			else if( option == "FD_ORDER_FOUR" )
				target().fd_method_order( target_t::FD_ORDER_FOUR );
			else if( option == "FD_ORDER_FOUR_CENTRAL" )
				target().fd_method_order( target_t::FD_ORDER_FOUR_CENTRAL );
			else if( option == "FD_ORDER_FOUR_AUTO" )
				target().fd_method_order( target_t::FD_ORDER_FOUR_AUTO );
			else
				TEST_FOR_EXCEPTION(
					true, std::invalid_argument
					,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
					"\"fd_method_order\".  Only the options FD_ORDER_ONE, FD_ORDER_TWO, "
					"FD_ORDER_TWO_CENTRAL, FD_ORDER_TWO_AUTO, FD_ORDER_FOUR, FD_ORDER_FOUR_CENTRAL "
					"and FD_ORDER_FOUR_AUTO are available" );
			break;
		}
	    case FD_STEP_SELECT:
		{
			const std::string &option = option_value.c_str();
			if( option == "FD_STEP_ABSOLUTE" )
				target().fd_step_select( target_t::FD_STEP_ABSOLUTE );
			else if( option == "FD_STEP_RELATIVE" )
				target().fd_step_select( target_t::FD_STEP_RELATIVE );
			else
				TEST_FOR_EXCEPTION(
					true, std::invalid_argument
					,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
					"\"fd_step_select\".  Only the options are available" );
			break;
		}
	    case FD_STEP_SIZE:
			target().fd_step_size(::atof(option_value.c_str()));
			break;
	    case FD_STEP_SIZE_MIN:
			target().fd_step_size_min(::atof(option_value.c_str()));
			break;
	    case FD_STEP_SIZE_F:
			target().fd_step_size_f(::atof(option_value.c_str()));
			break;
	    case FD_STEP_SIZE_C:
			target().fd_step_size_c(::atof(option_value.c_str()));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace NLPInterfacePack
