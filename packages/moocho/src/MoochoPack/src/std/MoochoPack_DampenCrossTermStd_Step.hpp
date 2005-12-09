// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_DampenCrossTermStd_Step.hpp
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

#ifndef DAMPEN_CROSS_TERM_STD_STEP_H
#define DAMPEN_CROSS_TERM_STD_STEP_H

#include "../rSQPAlgo_Step.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Compute a dampening term zeta_k for the cross term w_k such that
 * Gf'*Z*pz <= 0.
 *
 * This condition Gf'*Z*pz <= 0 is needed to ensure descent of many
 * merit functions.
 *
 * This implementation ensures that Gf'*Z*pz <= 0 only if
 * there will not be any active constraints when the reduced QP subproblem
 * is solved (nu_k = 0) and there are no undecomposed equality constraints
 * or if there they are linearly dependent (lambda_k(undecomp_con) = 0).
 * 
 * In particular this implementation computes zeta_k such that:
 * 
 * Gf'*Z*pz <= frac_descent * rGf'inv(B)*rGf
 * 
 * where: 0 < frac_descent < 1
 * 
 * To ensure strong descent (and hopefully deal with the cases where
 * nu_k != 0 and lambda_k(undecomp_con) != 0) the parameter frac_descent
 * is set to frac_descent = 0.9 by default.
 * 
 * The basis derivation goes like this:
 * 
 * ToDo: Finish documentation!
 */
class DampenCrossTermStd_Step : public rSQPAlgo_Step {
public:

	/// «std comp» members for frac_descent
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, frac_descent )

	///
	DampenCrossTermStd_Step(const value_type& frac_descent = 0.9);

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class DampenCrossTermStd_Step

}	// end namespace MoochoPack 

#endif	// DAMPEN_CROSS_TERM_STD_STEP_H
