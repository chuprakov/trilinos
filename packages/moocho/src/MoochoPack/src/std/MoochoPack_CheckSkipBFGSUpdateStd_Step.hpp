// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_CheckSkipBFGSUpdateStd_Step.hpp
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

#ifndef CHECK_SKIP_BFGS_UPDATE_STD_STEP_H
#define CHECK_SKIP_BFGS_UPDATE_STD_STEP_H

#include "MoochoPack_quasi_newton_stats.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Checks if a BFGS update should be preformed.
  */
class CheckSkipBFGSUpdateStd_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	///
	/** <<std member comp>> members for proportionality constant to use in the
	  * test for if to perform BFGS update.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, skip_bfgs_prop_const )

	///
	CheckSkipBFGSUpdateStd_Step(
		value_type	skip_bfgs_prop_const	= 10.0
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
	quasi_newton_stats_iq_member	quasi_newton_stats_;
	
};	// end class ReducedHessianBFGS_Step

}	// end namespace MoochoPack 

#endif	// CHECK_SKIP_BFGS_UPDATE_STD_STEP_H
