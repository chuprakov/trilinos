// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgo.h
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

#ifndef RSQP_ALGO_H
#define RSQP_ALGO_H

#include "rSQPAlgoInterface.h"
#include "rSQPAlgoContainer.h"
#include "rSQPState.h"
#include "GeneralIterationPack/src/Algorithm.h"
#include "StandardAggregationMacros.h"

namespace ReducedSpaceSQPPack {

///
/** rSQP Algorithm control class.
  */
class rSQPAlgo
	: public rSQPAlgoInterface
	, public GeneralIterationPack::Algorithm
{
public:

	/** @name Public Types */
	//@{

	//@}

	/// Constructs with no step, added_step, pre_step, post_step, state, or decomp_sys objects added.
	rSQPAlgo();

	/// <<std aggr>> members for algo_cntr
	STANDARD_AGGREGATION_MEMBERS( rSQPAlgoContainer, algo_cntr )

	/// <<std aggr>> members for nlp
	STANDARD_AGGREGATION_MEMBERS( NLP, nlp )

	///
	rSQPState& rsqp_state()
	{	return dynamic_cast<rSQPState&>(state()); }

	///
	const rSQPState& rsqp_state() const
	{	return dynamic_cast<const rSQPState&>(state()); }

	///
	void do_step_first(Algorithm::poss_type first_step_poss)
	{	first_step_poss_ = first_step_poss; }

	/** @name Overridden form rSQPAlgoInteface */
	//@{	
	
	///
	const rSQPState& retrieve_state() const;

	///
	/** This is the main control function for the rSQP algorithm.
	  *
	  * This function basically just calls Algorithm::do_algorithm(...).
	  */
	rSQPSolverClientInterface::EFindMinReturn dispatch();

	///
	void interface_print_algorithm(std::ostream& out) const;
	///
	void interface_set_algo_timing( bool algo_timing );
	///
	bool interface_algo_timing() const;
	///
	void interface_print_algorithm_times( std::ostream& out ) const;

	//@}

	/// overridden from Algorihth.

	///
	void print_algorithm(std::ostream& out) const;

protected:

	// First step to execute
	Algorithm::poss_type first_step_poss_;

};	// end class rSQPAlgo

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_H
