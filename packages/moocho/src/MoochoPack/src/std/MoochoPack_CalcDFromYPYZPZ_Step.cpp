// ////////////////////////////////////////////////////////////////////////////
// CalcDFromYPYZPZ_Step.cpp
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

#include <limits>
#include <ostream>

#include "ReducedSpaceSQPPack/src/std/CalcDFromYPYZPZ_Step.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "GeneralIterationPack/src/print_algorithm_step.hpp"
//#include "ConstrainedOptimizationPack/src/print_vector_change_stats.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::CalcDFromYPYZPZ_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using AbstractLinAlgPack::dot;
	using LinAlgOpPack::V_VpV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// d = Ypy + Zpz
	VectorWithOpMutable    &d_k    = s.d().set_k(0);
	const VectorWithOp     &Ypy_k = s.Ypy().get_k(0);
	const VectorWithOp     &Zpz_k = s.Zpz().get_k(0);
	V_VpV( &d_k, Ypy_k, Zpz_k );

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		const value_type very_small = std::numeric_limits<value_type>::min();
		out << "\n(Ypy_k'*Zpz_k)/(||Ypy_k||2 * ||Zpz_k||2 + eps) = "
			<< dot( Ypy_k, Zpz_k ) / ( Ypy_k.norm_2() * Zpz_k.norm_2() + very_small );
		out	<< "\n||d||inf = " << d_k.norm_inf() << std::endl;
/*
        ConstrainedOptimizationPack::print_vector_change_stats(
            s.x().get_k(0), "x", s.d().get_k(0), "d", out );
*/
		// ToDo: Replace the above with a reduction operator!
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nd_k = \n" << d_k;
	}

	return true;
}

void ReducedSpaceSQPPack::CalcDFromYPYZPZ_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculates the search direction d from Ypy and Zpz\n"
		<< L << "d_k = Ypy_k + Zpz_k \n";
}
