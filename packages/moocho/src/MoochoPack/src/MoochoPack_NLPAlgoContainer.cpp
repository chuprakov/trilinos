// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoContainer.cpp

//#define RELEASE_TRACE

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <assert.h>

#include <iostream>	// used for debugging the Release version.

#include "Misc/include/debug.h"
#include "../include/rSQPAlgoContainer.h"
#include "../include/rSQPAlgoInterface.h"

void ReducedSpaceSQPPack::rSQPAlgoContainer::set_max_run_time(double max_run_time)
{
	algo().set_max_run_time(max_run_time);
}

double ReducedSpaceSQPPack::rSQPAlgoContainer::max_run_time() const
{
	return algo().return_max_run_time();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::set_config(const config_ptr_t& config)
{
#ifdef RELEASE_TRACE
	std::cout << "\n*** if( get_config().get() ) this->config().deconfig_algo_cntr(*this) ...\n";
#endif
	if( get_config().get() ) this->config().deconfig_algo_cntr(*this);
	config_ = config;
#ifdef RELEASE_TRACE
	std::cout << "\n*** if( get_config().get() ) this->config().config_algo_cntr(*this) ...\n";
#endif
	if( get_config().get() ) this->config().config_algo_cntr(*this);	// The old way
}

ReducedSpaceSQPPack::rSQPSolverClientInterface::EFindMinReturn
ReducedSpaceSQPPack::rSQPAlgoContainer::find_min() {
//	TRACE0( "\n*** rSQPAlgoContainer::find_min() called\n" );

//	configure_algorithm();

//	TRACE0( "\n*** After configure_algorithm() was called\n" );
//	config().print_state();

	config().init_algo(algo());
	return algo().dispatch();
}

const ReducedSpaceSQPPack::rSQPState& ReducedSpaceSQPPack::rSQPAlgoContainer::state() const {
	return algo().retrieve_state();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::configure_algorithm() {
//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::configure_algorithm() called\n" );

	assert_valid_setup();
	if(!get_algo())
		config().config_algo_cntr(*this);

//	TRACE0( "\n*** After config().config_algo_cntr(*this) was called\n" );
//	config().print_state();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::print_algorithm(std::ostream& out) const {
//	TRACE0( "\n*** rSQPAlgoContainer::print_algorithm(out) called\n" );
//	TRACE0( "\n*** Before config().print_state() was called\n" );
//	config().print_state();
	assert(get_algo());
	algo().interface_print_algorithm(out);
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::set_algo_timing( bool algo_timing )
{
	assert(get_algo());
	algo().interface_set_algo_timing(algo_timing);
}

bool ReducedSpaceSQPPack::rSQPAlgoContainer::algo_timing() const
{
	assert(get_algo());
	return algo().interface_algo_timing();
}

void ReducedSpaceSQPPack::rSQPAlgoContainer::print_algorithm_times(
	std::ostream& out ) const
{
	assert(get_algo());
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