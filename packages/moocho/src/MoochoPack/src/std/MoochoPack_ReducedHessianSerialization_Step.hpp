// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianSerialization_Step.hpp
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

#ifndef REDUCED_HESSIAN_SERIALIZATION_STEP_HPP
#define REDUCED_HESSIAN_SERIALIZATION_STEP_HPP

#include "MoochoPack/src/std/ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack/src/std/quasi_newton_stats.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Serializes rHL_k to and from a file.
 *
 * ToDo: Finish documentation!
 */
class ReducedHessianSerialization_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/// Pick the file name to read in the reduced Hessian from
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, reduced_hessian_input_file_name );

	/// Pick the file name to write in the reduced Hessian to
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, reduced_hessian_output_file_name );

	///
	ReducedHessianSerialization_Step(
		const std::string    &reduced_hessian_input_file_name   = "reduced_hessian.in"
		,const std::string   &reduced_hessian_output_file_name  = "reduced_hessian.out"
		);
	
	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(
		Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		,poss_type assoc_step_poss
		);
	///
  void finalize_step(
		Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		,poss_type assoc_step_poss
		);
	///
	void print_step(
		const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		,poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str
		) const;
	//@}
	
};	// end class ReducedHessianSerialization_Step

}	// end namespace MoochoPack 

#endif	// REDUCED_HESSIAN_SERIALIZATION_STEP_HPP
