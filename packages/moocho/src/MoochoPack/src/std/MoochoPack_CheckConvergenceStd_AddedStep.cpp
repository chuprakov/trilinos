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

#include "MoochoPack/src/std/CheckConvergenceStd_AddedStep.hpp"
#include "MoochoPack/src/NLPAlgoContainer.hpp"
#include "MoochoPack/src/moocho_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"

namespace MoochoPack {

CheckConvergenceStd_AddedStep::CheckConvergenceStd_AddedStep(
  Teuchos::RefCountPtr<CheckConvergence_Strategy> convergence_strategy
	)
	:
	convergence_strategy_(convergence_strategy)
	{}

bool CheckConvergenceStd_AddedStep::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
	{
	
    TEST_FOR_EXCEPTION(!convergence_strategy_.get(),
					std::logic_error,
					"Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
	  );

	NLPAlgo	&algo	  = rsqp_algo(_algo);

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	const bool found_solution = convergence_strategy_->Converged(_algo);

	if( found_solution )
	{
		if( static_cast<int>(olevel) > static_cast<int>(PRINT_NOTHING) )
			out	<< "\nJackpot!  Found the solution!!!!!! (k = " << algo.state().k() << ")\n";
		algo.terminate(true);	// found min
		return false; // skip the other steps and terminate
	}
	
	if( static_cast<int>(olevel) > static_cast<int>(PRINT_NOTHING) )
		out	<< "\nHave not found the solution yet, have to keep going :-(\n";
	
	// We are not at the solution so keep going
	return true;
	}

void CheckConvergenceStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
	{

    TEST_FOR_EXCEPTION(!convergence_strategy_.get(),
					std::logic_error,
					"Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
	  );
	
	convergence_strategy_->print_step(algo, out, L);
	}

}	// end namespace MoochoPack
