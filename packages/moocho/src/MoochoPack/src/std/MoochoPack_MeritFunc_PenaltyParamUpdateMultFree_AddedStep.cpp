// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateMultFree_AddedStep.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <math.h>

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdateMultFree_AddedStep::MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
	value_type    small_mu
	,value_type   mult_factor
	,value_type   kkt_near_sol
	)
	:MeritFunc_PenaltyParamUpdateGuts_AddedStep(small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateMultFree_AddedStep::min_mu(
	rSQPState& s, value_type* min_mu
	) const
{
	using AbstractLinAlgPack::dot;

	IterQuantityAccess<VectorWithOpMutable>
		&Gf_iq    = s.Gf(),
		&nu_iq    = s.nu(),
		&Ypy_iq   = s.Ypy(),
		&c_iq     = s.c();
	if ( Gf_iq.updated_k(0) && nu_iq.updated_k(0) && Ypy_iq.updated_k(0) && c_iq.updated_k(0) ) {
		// min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,1)
		*min_mu = ::fabs( 	  dot( Gf_iq.get_k(0), Ypy_iq.get_k(0) )
							+ dot( nu_iq.get_k(0), Ypy_iq.get_k(0) )
						) / c_iq.get_k(0).norm_1();
		return true;
	}
	return false;
}

void MeritFunc_PenaltyParamUpdateMultFree_AddedStep::print_min_mu_step(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if Gf_k, nu_k, Ypy_k and c_k are updated then\n"
		<< L << "   min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,1)\n"
		<< L << "   update_mu = true\n"
		<< L << "else\n"
		<< L << "   update_mu = false\n"
		<< L << "endif\n"
		;
}

}	// end namespace ReducedSpaceSQPPack
