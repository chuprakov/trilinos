// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateStd_Step.hpp
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

#include "MoochoPack/src/std/ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack/src/std/quasi_newton_stats.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Updates rHL_k using a secant update.
 *
 * This class will initialize rHL = I if it is not already and will handle
 * transitions to new basis selections by reseting rHL = I.  The actually
 * update is delegated to a strategy object (see #secant_update# below).
 *
 * See the printed step documentation for a description.
 */
class ReducedHessianSecantUpdateStd_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	///
	/** <<std comp>> members for the strategy interface object that will
	 * actually perform the secant update.
	 */
	STANDARD_COMPOSITION_MEMBERS( ReducedHessianSecantUpdate_Strategy, secant_update )

	///
	ReducedHessianSecantUpdateStd_Step(
		const secant_update_ptr_t& secant_update = Teuchos::null
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
	enum { NO_BASIS_UPDATED_YET = INT_MIN };
	int                             num_basis_;
	int                             iter_k_rHL_init_ident_;
	quasi_newton_stats_iq_member	quasi_newton_stats_;
	
};	// end class ReducedHessianSecantUpdateStd_Step

}	// end namespace MoochoPack 

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H
