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

#include "MoochoPack_MeritFunc_DummyUpdate_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

bool MeritFunc_DummyUpdate_Step::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	)
{
	using Teuchos::dyn_cast;

	NLPAlgo	&algo	= rsqp_algo(_algo);
	NLPAlgoState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
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
	const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Simply sets the current merit function value for unconstrained linesearch\n"
		<< L << "if merit_func_nlp_k not updated set merit_func_nlp_k = merit_func_nlp_k(last_updated)\n";
}

} // end namespace MoochoPack
