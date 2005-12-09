// ////////////////////////////////////////////////////////////////////////////
// NumFixedDepIndep_AddedStep.cpp
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

#include "../std/MoochoPack_NumFixedDepIndep_AddedStep.hpp"
#include "../MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"

bool MoochoPack::NumFixedDepIndep_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	NLPAlgo	&algo	= rsqp_algo(_algo);
	NLPAlgoState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( s.nu().updated_k(0) && s.nu().get_k(0).nz() ) {
		const Range1D
			dep		= s.var_dep(),
			indep	= s.var_indep();
		const SpVector &nu_k	= s.nu().get_k(0);
		size_type fixed_dep = 0, fixed_indep = 0;
		for( SpVector::const_iterator itr = nu_k.begin(); itr != nu_k.end(); ++itr ) {
			if( dep.in_range( itr->indice() + nu_k.offset() ) )
				fixed_dep++;
			else if( indep.in_range( itr->indice() + nu_k.offset() ) )
				fixed_indep++;
			else
				assert(0);	// should never happen
		}
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nnum_dep_fixed = "		<< fixed_dep
				<< "\nnum_indep_fixed = "	<< fixed_indep << std::endl;
		}
	}
	else {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nnu not calculated for the kth iteration\n";
		}
	}

	return true;
}

void MoochoPack::NumFixedDepIndep_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Counts the number of fixed variables from "
				"the dependent and independent sets\n";
}
