// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgo.cpp

// #define RELEASE_TRACE

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <iostream> // for debugging Release version.
#include <typeinfo>

#include "../include/rSQPAlgo.h"

namespace ReducedSpaceSQPPack {

void rSQPAlgo::set_algo_cntr(rSQPAlgoContainer* algo_cntr) {
#ifdef RELEASE_TRACE
	std::cout << "\n*** algo_cntr_ = ";
	std::cout << algo_cntr_ << "\n";	// Accessing this causes an unknown exception to be thrown.
	std::cout << "\n*** before algo_cntr_ = algo_cntr ...\n";
#endif
	algo_cntr_ = algo_cntr;
#ifdef RELEASE_TRACE
	std::cout << "\n*** after algo_cntr_ = algo_cntr ...\n";
#endif
}

// Overridden form rSQPAlgoInteface

const rSQPState& rSQPAlgo::retrieve_state() const
{
	return dynamic_cast<const rSQPState&>(state());
}

void rSQPAlgo::set_max_run_time(double max_run_time)
{
	this->max_run_time(max_run_time);
}

double rSQPAlgo::return_max_run_time() const
{
	return Algorithm::max_run_time();
}

rSQPSolverClientInterface::EFindMinReturn
rSQPAlgo::dispatch() {
	switch( do_algorithm(first_step_poss_) ) {
		case GeneralIterationPack::TERMINATE_TRUE:
			return rSQPSolverClientInterface::SOLUTION_FOUND;
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
		<< "\n*** NLPReduced ***\n"
		<< typeid(*get_nlp()).name() << "\n";

	Algorithm::print_algorithm(out);
}

}	// end namespace ReducedSpaceSQPPack
