// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateWithMult_AddedStep.cpp
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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateWithMult_AddedStep::MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
		  const merit_func_ptr_t& merit_func, value_type small_mu
		, value_type mult_factor, value_type kkt_near_sol )
	: MeritFunc_PenaltyParamUpdateGuts_AddedStep(merit_func,small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateWithMult_AddedStep::min_mu(
	NLPAlgoState& s, value_type* min_mu ) const
{
	if ( s.lambda().updated_k(0) ) {
		*min_mu = s.lambda().get_k(0).norm_inf();
		return true;
	}
	return false;
}

void MeritFunc_PenaltyParamUpdateWithMult_AddedStep::print_min_mu_step(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if lambda_k is updated then\n"
		<< L << "    min_mu = norm( lambda_k, inf )\n"
		<< L << "    update_mu = true\n"
		<< L << "else\n"
		<< L << "    update_mu = false\n"
		<< L << "endif\n"
		;
}

}	// end namespace MoochoPack
