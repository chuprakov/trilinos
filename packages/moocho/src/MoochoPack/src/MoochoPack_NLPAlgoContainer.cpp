// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoContainer.cpp

//#define RELEASE_TRACE

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <assert.h>

#include <iostream>	// used for debugging the Release version.

#include "../include/rSQPAlgoContainer.h"
#include "../include/rSQPAlgoInterface.h"
#include "Misc/include/debug.h"

void ReducedSpaceSQPPack::rSQPAlgoContainer::max_run_time(double max_run_time)
{
	algo().set_max_run_time(max_run_time);
}

double ReducedSpaceSQPPack::rSQPAlgoContainer::max_run_time() const
{
	return algo().return_max_run_time();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::set_config(const config_ptr_t& config)
{
	algo_ = algo_ptr_t(0); // Remove our reference to the current (configured?) algorithm.
	config_ = config;
}

ReducedSpaceSQPPack::rSQPSolverClientInterface::EFindMinReturn
ReducedSpaceSQPPack::rSQPAlgoContainer::find_min() {
	config().init_algo(algo());
	return algo().dispatch();
}

const ReducedSpaceSQPPack::rSQPState& ReducedSpaceSQPPack::rSQPAlgoContainer::state() const {
	return algo().retrieve_state();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::configure_algorithm(std::ostream* trase_out) {
	assert_valid_setup();
	if(!get_algo().get())
		config().config_algo_cntr(*this,trase_out);
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::print_algorithm(std::ostream& out) const {
	algo().interface_print_algorithm(out);
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::set_algo_timing( bool algo_timing )
{
	algo().interface_set_algo_timing(algo_timing);
}

bool ReducedSpaceSQPPack::rSQPAlgoContainer::algo_timing() const
{
	return algo().interface_algo_timing();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::print_algorithm_times(
	std::ostream& out ) const
{
	algo().interface_print_algorithm_times(out);
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::assert_valid_setup() const {
	if( !get_nlp().get() )
		throw rSQPSolverClientInterface::InvalidSetup( "rSQPAlgoContainer::assert_valid_setup() : "
				"The NLPReduced object has not been set" );
	if( !get_track().get() )
		throw rSQPSolverClientInterface::InvalidSetup( "rSQPAlgoContainer::assert_valid_setup() : "
				"The rSQPTrack object has not been set" );
	if( !get_config().get() )
		throw rSQPSolverClientInterface::InvalidSetup( "rSQPAlgoContainer::assert_valid_setup() : "
				"The rSQPAlgo_Config object has not been set" );
}