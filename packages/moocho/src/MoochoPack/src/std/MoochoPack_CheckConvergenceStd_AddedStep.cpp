// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.cpp
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

//#include <assert.h>

#include <ostream>

#include "ReducedSpaceSQPPack/src/std/CheckConvergenceStd_AddedStep.h"
#include "ReducedSpaceSQPPack/src/rSQPAlgoContainer.h"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"

namespace ReducedSpaceSQPPack {

CheckConvergenceStd_AddedStep::CheckConvergenceStd_AddedStep(
  MemMngPack::ref_count_ptr<CheckConvergence_Strategy> convergence_strategy
	)
	:
	convergence_strategy_(convergence_strategy)
	{}

bool CheckConvergenceStd_AddedStep::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
	{
	
    THROW_EXCEPTION(!convergence_strategy_.get(),
					std::logic_error,
					"Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
	  );

	rSQPAlgo	&algo	  = rsqp_algo(_algo);
	rSQPState	&s		  = algo.rsqp_state();
	NLP			&nlp	  = algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	const bool found_solution = convergence_strategy_->Converged(_algo);

	if( found_solution )
		{
		const size_type
			m  = nlp.m(),
			mI = nlp.mI();
		
		IterQuantityAccess<VectorWithOpMutable>
			&x_iq       = s.x(),
			*lambda_iq  = m  ? &s.lambda()  : NULL,
			*lambdaI_iq = mI ? &s.lambdaI() : NULL,
			&nu_iq      = s.nu();
		
		nlp.report_final_solution(
		  x_iq.get_k(0),
		  m  && lambda_iq->updated_k(0)  ? &lambda_iq->get_k(0)  : NULL,
		  mI && lambdaI_iq->updated_k(0) ? &lambdaI_iq->get_k(0) : NULL,
		  nu_iq.updated_k(0)            ? &nu_iq.get_k(0)        : NULL,
		  true
		  );

		out	<< "\nJackpot!  Found the solution!!!!!! (k = " << algo.state().k() << ")\n";
		algo.terminate(true);	// found min
		return false; // skip the other steps and terminate
		}

	out	<< "\nHave not found the solution yet, have to keep going :-(\n";
	
	// We are not at the solution so keep going
	return true;
	}

void CheckConvergenceStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
	{

    THROW_EXCEPTION(!convergence_strategy_.get(),
					std::logic_error,
					"Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
	  );
	
	convergence_strategy_->print_step(algo, out, L);
	}

}	// end namespace ReducedSpaceSQPPack
