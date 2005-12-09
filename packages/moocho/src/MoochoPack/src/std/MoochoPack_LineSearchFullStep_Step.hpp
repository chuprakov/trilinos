// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_LineSearchFullStep_Step.hpp
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

#ifndef LINE_SEARCH_FULL_STEP_STEP_H
#define LINE_SEARCH_FULL_STEP_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "ConstrainedOptPack_VariableBoundsTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Takes the full step <tt>x_kp1 = x_k + d_k (d_k = Ypy_k + Zpz_k)</tt>.
  */
class LineSearchFullStep_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/// «std comp» Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )

	///
	LineSearchFullStep_Step(
		const bounds_tester_ptr_t&	bounds_tester
		);

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:
	/// Not defined and not to be called
	LineSearchFullStep_Step();

};	// end class LineSearchFullStep_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FULL_STEP_STEP_H
