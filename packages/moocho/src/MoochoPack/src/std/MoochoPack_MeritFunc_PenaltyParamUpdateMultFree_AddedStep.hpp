// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp
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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.hpp"

namespace MoochoPack {

///
/** Specializes the update of the penalty parameter for a merit function as:
  * min_mu = |(Gf_k+nu_k)'* Ypy_k| / ||c_k||1.
  */
class MeritFunc_PenaltyParamUpdateMultFree_AddedStep
	: public MeritFunc_PenaltyParamUpdateGuts_AddedStep
{
public:

	///
	MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
		value_type   small_mu = 1e-6
		,value_type  mult_factor = 1e-4
		,value_type  kkt_near_sol = 1.0
		);

protected:

	/** @name Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep */
	//@{

	///
	bool min_mu( NLPAlgoState& s, value_type* min_mu ) const;
	///
	void print_min_mu_step( std::ostream& out
		, const std::string& leading_str ) const;

	//@}
	
};	// end class MeritFunc_PenaltyParamUpdateMultFree_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
