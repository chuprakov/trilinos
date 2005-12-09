// ////////////////////////////////////////////////////////////////////////////
// NLPAlgo.cpp
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

// #define RELEASE_TRACE

#include <iostream> // for debugging Release version.
#include <typeinfo>

#include "MoochoPack_NLPAlgo.hpp"

namespace MoochoPack {

NLPAlgo::NLPAlgo()
	: algo_cntr_(NULL), nlp_(NULL), first_step_poss_(1)
{}

// Overridden form rSQPAlgoInteface

const NLPAlgoState& NLPAlgo::retrieve_state() const
{
	return dynamic_cast<const NLPAlgoState&>(state());
}

NLPSolverClientInterface::EFindMinReturn
NLPAlgo::dispatch() {
	switch( do_algorithm(first_step_poss_) ) {
		case IterationPack::TERMINATE_TRUE:
			return NLPSolverClientInterface::SOLUTION_FOUND;
		case IterationPack::TERMINATE_FALSE:
			return NLPSolverClientInterface::ALGORITHMIC_ERROR;
		case IterationPack::MAX_ITER_EXCEEDED:
			return NLPSolverClientInterface::MAX_ITER_EXCEEDED;
		case IterationPack::MAX_RUN_TIME_EXCEEDED:
			return NLPSolverClientInterface::MAX_RUN_TIME_EXCEEDED;
		case IterationPack::INTERRUPTED_TERMINATE_TRUE:
			return NLPSolverClientInterface::SOLUTION_FOUND;
		case IterationPack::INTERRUPTED_TERMINATE_FALSE:
			return NLPSolverClientInterface::ALGORITHMIC_ERROR;
		default:
			assert(0);
	}
	return NLPSolverClientInterface::SOLUTION_FOUND;	// will never be called.
}

void NLPAlgo::interface_print_algorithm(std::ostream& out) const {
	print_steps(out);
	print_algorithm(out);
}

void NLPAlgo::interface_set_algo_timing( bool algo_timing ) {
	set_algo_timing(algo_timing);
}

bool NLPAlgo::interface_algo_timing() const {
	return algo_timing();
}

void NLPAlgo::interface_print_algorithm_times( std::ostream& out ) const {
	print_algorithm_times(out);
}

// Overridden from Algorithm.

void NLPAlgo::print_algorithm(std::ostream& out) const {
	out
		<< "\n*** NLP ***\n"
		<< typeid(*get_nlp()).name() << "\n";

	Algorithm::print_algorithm(out);
}

}	// end namespace MoochoPack
