// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoContainer.cpp
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

//#define RELEASE_TRACE

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <assert.h>

#include <iostream>	// used for debugging the Release version.

#include "ReducedSpaceSQPPack/include/rSQPAlgoContainer.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgoInterface.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "NLPInterfacePack/include/NLP.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "ThrowException.h"

namespace {

void report_final_failure( const ReducedSpaceSQPPack::rSQPState& s, NLPInterfacePack::NLP* nlp )
{
	const AbstractLinAlgPack::size_type
		m  = nlp->m(),
		mI = nlp->mI(),
		nb = nlp->num_bounded_x();
	nlp->report_final_solution(
		s.x().get_k(0)                                                     // x
		,( m  && s.lambda().updated_k(0) ) ? &s.lambda().get_k(0)  : NULL  // lambda
		,( mI && s.lambda().updated_k(0) ) ? &s.lambdaI().get_k(0) : NULL  // lambdaI
		,( nb && s.nu().updated_k(0)     ) ? &s.nu().get_k(0)      : NULL  // nu
		, false                                                            // optimal = false
		);
}

} // end namespace

namespace ReducedSpaceSQPPack {

// Overridden from rSQPAlgoClient interface

void rSQPAlgoContainer::set_config(const config_ptr_t& config)
{
	algo_ = algo_ptr_t(NULL); // Remove our reference to the current (configured?) algorithm.
	config_ = config;
}

rSQPAlgoContainer::config_ptr_t&
rSQPAlgoContainer::get_config()
{	
	return config_;
}

const rSQPAlgoContainer::config_ptr_t&
rSQPAlgoContainer::get_config() const
{	
	return config_;
}

rSQPAlgo_Config&
rSQPAlgoContainer::config()
{	
	return *config_;
}

const rSQPAlgo_Config&
rSQPAlgoContainer::config() const
{	
	return *config_;
}

rSQPSolverClientInterface::EFindMinReturn
rSQPAlgoContainer::find_min()
{
	config().init_algo(&algo());
	EFindMinReturn solve_return;
	try {
		solve_return = algo().dispatch();
	}
	catch(...) {
		report_final_failure(algo().retrieve_state(),&nlp());
		throw;
	}
	if( solve_return != SOLUTION_FOUND ) {
		report_final_failure(algo().retrieve_state(),&nlp());
	}
	return solve_return;
}

void rSQPAlgoContainer::configure_algorithm(std::ostream* trase_out)
{
	assert_valid_setup();
	if(!get_algo().get())
		config().config_algo_cntr(this,trase_out);
}

void rSQPAlgoContainer::print_algorithm(std::ostream& out) const
{
	algo().interface_print_algorithm(out);
}

void rSQPAlgoContainer::set_algo_timing( bool algo_timing )
{
	algo().interface_set_algo_timing(algo_timing);
}

bool rSQPAlgoContainer::algo_timing() const
{
	return algo().interface_algo_timing();
}

void rSQPAlgoContainer::print_algorithm_times(
	std::ostream& out ) const
{
	algo().interface_print_algorithm_times(out);
}

void rSQPAlgoContainer::assert_valid_setup() const {
	THROW_EXCEPTION(
		!get_nlp().get(), rSQPSolverClientInterface::InvalidSetup
		,"rSQPAlgoContainer::assert_valid_setup() : The NLP object has not been set" );
	THROW_EXCEPTION(
		!get_track().get(), rSQPSolverClientInterface::InvalidSetup
		,"rSQPAlgoContainer::assert_valid_setup() : The rSQPTrack object has not been set" );
	THROW_EXCEPTION(
		!get_config().get(), rSQPSolverClientInterface::InvalidSetup
		,"rSQPAlgoContainer::assert_valid_setup() : The rSQPAlgo_Config object has not been set" );
}

} // end namespace ReducedSpaceSQPPack
