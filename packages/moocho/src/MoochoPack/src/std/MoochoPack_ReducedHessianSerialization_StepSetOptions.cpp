// ////////////////////////////////////////////////////////////////
// ReducedHessianSerialization_StepSetOptions.cpp
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

#include "MoochoPack_ReducedHessianSerialization_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

const int local_num_options = 2;

enum local_EOptions {
	REDUCED_HESSIAN_INPUT_FILE_NAME
	,REDUCED_HESSIAN_OUTPUT_FILE_NAME
};

const char* local_SOptions[local_num_options]	= {
	"reduced_hessian_input_file_name"
	,"reduced_hessian_output_file_name"
};

inline
std::string remove_quotes( const std::string &option_name, const std::string &str )
{
	TEST_FOR_EXCEPTION(
		str[0]!='\"' || str[str.length()-1]!='\"', std::logic_error
		,"Error, the option \'" << option_name << "\' must have a single set of quotes around it!"
		);
	if(str.length()==2)
		return std::string("");
	return str.substr(1,str.length()-2);
}

} // namespace

namespace MoochoPack {

ReducedHessianSerialization_StepSetOptions::ReducedHessianSerialization_StepSetOptions(
			  ReducedHessianSerialization_Step* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			ReducedHessianSerialization_Step >( target )
{}

void ReducedHessianSerialization_StepSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;

	typedef ReducedHessianSerialization_Step target_t;
	switch( (local_EOptions)option_num ) {
		case REDUCED_HESSIAN_INPUT_FILE_NAME : {
			target().reduced_hessian_input_file_name(remove_quotes("reduced_hessian_input_file_name",option_value));
			break;
		}
		case REDUCED_HESSIAN_OUTPUT_FILE_NAME : {
			target().reduced_hessian_output_file_name(remove_quotes("reduced_hessian_output_file_name",option_value));
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
