// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgo.cpp
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

#include "../include/rSQPAlgo.h"

namespace ReducedSpaceSQPPack {

rSQPAlgo::rSQPAlgo()
	: algo_cntr_(NULL), nlp_(NULL), first_step_poss_(1)
{}

// Overridden form rSQPAlgoInteface

const rSQPState& rSQPAlgo::retrieve_state() const
{
	return dynamic_cast<const rSQPState&>(state());
}

rSQPSolverClientInterface::EFindMinReturn
rSQPAlgo::dispatch() {
	switch( do_algorithm(first_step_poss_) ) {
		case GeneralIterationPack::TERMINATE_TRUE:
			return rSQPSolverClientInterface::SOLUTION_FOUND;
		case GeneralIterationPack::TERMINATE_FALSE:
			return rSQPSolverClientInterface::ALGORITHMIC_ERROR;
		case GeneralIterationPack::MAX_ITER_EXCEEDED:
			return rSQPSolverClientInterface::MAX_ITER_EXCEEDED;
		case GeneralIterationPack::MAX_RUN_TIME_EXCEEDED:
			return rSQPSolverClientInterface::MAX_RUN_TIME_EXCEEDED;
		default:
			assert(0);
	}
	return rSQPSolverClientInterface::SOLUTION_FOUND;	// will never be called.
}

void rSQPAlgo::interface_print_algorithm(std::ostream& out) const {
	print_steps(out);
	print_algorithm(out);
}

void rSQPAlgo::interface_set_algo_timing( bool algo_timing ) {
	set_algo_timing(algo_timing);
}

bool rSQPAlgo::interface_algo_timing() const {
	return algo_timing();
}

void rSQPAlgo::interface_print_algorithm_times( std::ostream& out ) const {
	print_algorithm_times(out);
}

// Overridden from Algorithm.

void rSQPAlgo::print_algorithm(std::ostream& out) const {
	out
		<< "\n*** NLP ***\n"
		<< typeid(*get_nlp()).name() << "\n";

	Algorithm::print_algorithm(out);
}

}	// end namespace ReducedSpaceSQPPack
