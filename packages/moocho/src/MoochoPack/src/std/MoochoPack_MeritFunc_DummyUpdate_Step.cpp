// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_DummyUpdate_Step.cpp
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

#include <limits>
#include <ostream>

#include "ReducedSpaceSQPPack/include/std/MeritFunc_DummyUpdate_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPDirecDeriv.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

bool MeritFunc_DummyUpdate_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	)
{
	using DynamicCastHelperPack::dyn_cast;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	IterQuantityAccess<MeritFuncNLP>
		&merit_func_nlp_iq = s.merit_func_nlp();

	if(!merit_func_nlp_iq.updated_k(0)) {
		const int last_updated_k = merit_func_nlp_iq.last_updated();
		MeritFuncNLP
			&merit_func_nlp_k =  ( last_updated_k != IterQuantity::NONE_UPDATED
								   ? merit_func_nlp_iq.set_k(0,last_updated_k)
								   : merit_func_nlp_iq.set_k(0) );
		MeritFuncNLPDirecDeriv
			&direc_deriv = dyn_cast<MeritFuncNLPDirecDeriv>(merit_func_nlp_k);

		direc_deriv.calc_deriv(
			s.Gf().get_k(0)
			,NULL
			,NULL
			,NULL
			,NULL
			,s.d().get_k(0)
			);

	}

	return true;
}

void MeritFunc_DummyUpdate_Step::print_step(
	const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Simply sets the current merit function value for unconstrained linesearch\n"
		<< L << "if merit_func_nlp_k not updated set merit_func_nlp_k = merit_func_nlp_k(last_updated)\n";
}

} // end namespace ReducedSpaceSQPPack
