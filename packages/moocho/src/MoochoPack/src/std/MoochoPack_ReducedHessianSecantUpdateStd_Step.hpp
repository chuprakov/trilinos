// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateStd_Step.h
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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H
#define REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Updates rHL_k using a secant update.
 *
 * This class will initialize rHL = I if it is not already and will handle
 * transitions to new basis selections by reseting rHL = I.  The actually
 * update is delegated to a strategy object (see #secant_update# below).
 *
 * See the printed step documentation for a description.
 */
class ReducedHessianSecantUpdateStd_Step : public ReducedHessian_Step {
public:

	///
	/** <<std comp>> members for the strategy interface object that will
	 * actually perform the secant update.
	 */
	STANDARD_COMPOSITION_MEMBERS( ReducedHessianSecantUpdate_Strategy, secant_update )

	///
	ReducedHessianSecantUpdateStd_Step(
		  const secant_update_ptr_t& secant_update = NULL
		);

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	enum { NO_BASIS_UPDATED_YET = -1 };
	int num_basis_;
	int iter_k_rHL_init_ident_;
	quasi_newton_stats_iq_member	quasi_newton_stats_;
	
};	// end class ReducedHessianSecantUpdateStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H
