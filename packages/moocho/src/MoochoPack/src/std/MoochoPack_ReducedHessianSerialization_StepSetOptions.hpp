// ////////////////////////////////////////////////////////////////
// MoochoPack_ReducedHessianSerialization_StepSetOptions.hpp
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

#ifndef REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H
#define REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSerialization_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for ReducedHessianSerialization_Step from an \c OptionsFromStream object.
 *
 * The default options group name is EvalNewPointStd.
 *
 * The options group is:
 *
 \verbatim
    options_group ReducedHessianSerialization {
		  reduced_hessian_input_file_name   = "reduced_hessian.in";
		  reduced_hessian_output_file_name  = "reduced_hessian.out";
    }
 \verbatim
 *
 * See the file Moocho.opt.NLPAlgoConfigMamaJama
 */
class ReducedHessianSerialization_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			ReducedHessianSerialization_Step >
{
public:

	///
	ReducedHessianSerialization_StepSetOptions(
		  ReducedHessianSerialization_Step* target = 0
		, const char opt_grp_name[] = "ReducedHessianSerialization" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSerialization_StepSetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H
