// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.hpp
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

#ifndef MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H

#include "MeritFunc_PenaltyParamUpdate_AddedStep.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncNLP.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Updates a set of penalty parameters for a merit function as:
  * mu(j) = max( mu(j), |lambda_k(j)| ).
  *
  * This class assumes the merit function supports the interfaces
  * MeritFuncPenaltyParams and MeritFuncNLPDirecDeriv.
  */
class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep
	: public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(const merit_func_ptr_t& merit_func
		, value_type small_mu = 1e-6, value_type min_mu_ratio = 1e-8
		, value_type mult_factor = 1e-4 , value_type kkt_near_sol = 1e-1 );

	// ///////////////////////////////
	// Overridden from AlgorithmStep

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss
		, IterationPack::EDoStepType type, poss_type assoc_step_poss
		, std::ostream& out, const std::string& leading_str ) const;

	// //////////////////////////////////////////////////////////////////////
	// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

	///
	void small_mu( value_type small_mu );
	///
	value_type small_mu() const;

	///
	void min_mu_ratio( value_type min_mu_ratio );
	///
	value_type min_mu_ratio() const;

	///
	void mult_factor( value_type mult_factor );
	///
	value_type mult_factor() const;

	///
	void kkt_near_sol( value_type kkt_near_sol );
	///
	value_type kkt_near_sol() const;

private:
	bool near_solution_;
	value_type small_mu_;
	value_type min_mu_ratio_;
	value_type mult_factor_;
	value_type kkt_near_sol_;
	value_type norm_inf_mu_last_;
	
};	// end class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
