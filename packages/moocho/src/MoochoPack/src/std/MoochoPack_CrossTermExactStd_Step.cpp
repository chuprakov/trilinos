// ////////////////////////////////////////////////////////////////////////////
// CrossTermExactStd_Step.cpp
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

#include "../std/CrossTermExactStd_Step.h"
#include "../rsqp_algo_conversion.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/src/VectorWithNorms.h"
#include "SparseLinAlgPack/src/MatrixWithOp.h"
#include "LinAlgPack/src/LinAlgOpPack.h"
#include "LinAlgPack/src/VectorClass.h"
#include "LinAlgPack/src/VectorOut.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::CrossTermExactStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgOpPack::V_MtV;
	using LinAlgPack::norm_inf;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	// tmp = HL * Ypy
	Vector tmp;
	V_MtV( &tmp, s.HL().get_k(0), BLAS_Cpp::no_trans, s.Ypy().get_k(0)() );
	// w = Z' * tmp
	V_MtV( &s.w().set_k(0).v(), s.Z().get_k(0), BLAS_Cpp::trans, tmp() );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||w||inf = "	<< s.w().get_k(0).norm_inf() << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nw_k =\n" << s.w().get_k(0)();
	}

	return true;
}

void ReducedSpaceSQPPack::CrossTermExactStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the exact reduced QP cross term\n"
		<< L << "w_k = Z_k' * HL_k * Ypy_k\n";
}
