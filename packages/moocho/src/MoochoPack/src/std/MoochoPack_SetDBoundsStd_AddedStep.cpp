// ////////////////////////////////////////////////////////////////////////////
// SetDBoundsStd_AddedStep.cpp
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

#include <ostream>

#include "ReducedSpaceSQPPack/include/std/SetDBoundsStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOut.h"

bool ReducedSpaceSQPPack::SetDBoundsStd_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Set the bounds on d.

	// dl = xl - x_k

	const Vector
		&x_k	= s.x().get_k(0).v();
	const SpVectorSlice
		&xl = algo.nlp().xl(),
		&xu = algo.nlp().xu();
	SparseBounds
		&d_bounds = d_bounds_(s).set_k(0);

	SpVector &dl = d_bounds.l;
	dl.uninitialized_resize( xl.size(), xl.nz(), xl.nz() );
	dl.assume_sorted(true);
	if(xl.nz()) {
		SpVectorSlice::const_iterator
			xl_itr		= xl.begin(),
			xl_itr_end	= xl.end();
		SpVector::iterator
			dl_itr		= dl.begin();
		for(; xl_itr != xl_itr_end; ++xl_itr, ++dl_itr) {
			const SpVectorSlice::element_type::indice_type i =  xl_itr->indice() + xl.offset();
			dl_itr->initialize( i, xl_itr->value() - x_k(i) );
		}
	}

	// du = xu - x_k

	SpVector &du = d_bounds.u;
	du.uninitialized_resize( xu.size(), xu.nz(), xu.nz() );
	du.assume_sorted(true);
	if(xu.nz()) {
		SpVectorSlice::const_iterator
			xu_itr		= xu.begin(),
			xu_itr_end	= xu.end();
		SpVector::iterator
			du_itr		= du.begin();
		for(; xu_itr != xu_itr_end; ++xu_itr, ++du_itr) {
			const SpVectorSlice::element_type::indice_type i =  xu_itr->indice() + xu.offset();
			du_itr->initialize( i, xu_itr->value() - x_k(i) );
		}
	}

	// Print out the QP bounds for the constraints
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nd_bounds_k.l = \n" << dl();
		out << "\nd_bounds_k.u = \n" << du();
	}

	if(algo.algo_cntr().check_results()) {
		dl.assert_valid_and_sorted();
		du.assert_valid_and_sorted();
	}

	return true;
}

void ReducedSpaceSQPPack::SetDBoundsStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Set the bounds on d\n"
		<< L << "d_bounds_k.l = xl - x_k\n"
		<< L << "d_bounds_k.u = xu - x_k\n"
		<< L << "if check_results == true then\n"
		<< L << "    assert that d_bounds_k.l and d_bounds_k.u are valid and sorted\n"
		<< L << "end\n"
		;
}
