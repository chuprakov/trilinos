// ////////////////////////////////////////////////////////////////////////////
// DampenCrossTermStd_Step.h

#ifndef DAMPEN_CROSS_TERM_STD_STEP_H
#define DAMPEN_CROSS_TERM_STD_STEP_H

#include "../rSQPAlgo_Step.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

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
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class DampenCrossTermStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// DAMPEN_CROSS_TERM_STD_STEP_H