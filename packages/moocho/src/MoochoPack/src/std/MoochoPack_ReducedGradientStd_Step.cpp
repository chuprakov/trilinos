// ////////////////////////////////////////////////////////////////////////////
// ReducedGradientStd_Step.cpp
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

#include "MoochoPack_ReducedGradientStd_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace MoochoPack {

bool ReducedGradientStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
	, poss_type assoc_step_poss
	)
{
	using BLAS_Cpp::trans;
	using LinAlgOpPack::V_MtV;

	NLPAlgo    &algo   = rsqp_algo(_algo);
	NLPAlgoState   &s      = algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	// Get iteration quantities
	IterQuantityAccess<VectorMutable>
		&Gf_iq  = s.Gf(),
		&rGf_iq = s.rGf();
	IterQuantityAccess<MatrixOp>
		&Z_iq = s.Z();

	// rGf = Z' * Gf
	V_MtV( &rGf_iq.set_k(0), Z_iq.get_k(0), trans, Gf_iq.get_k(0) );

	if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||rGf||inf = "	<< rGf_iq.get_k(0).norm_inf() << std::endl;
	}

	if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nrGf_k =\n" << rGf_iq.get_k(0);
	}

	return true;
}

void ReducedGradientStd_Step::print_step(
	const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Evaluate the reduced gradient of the objective funciton\n"
		<< L << "rGf_k = Z_k' * Gf_k\n";
}

} // end namespace MoochoPack
