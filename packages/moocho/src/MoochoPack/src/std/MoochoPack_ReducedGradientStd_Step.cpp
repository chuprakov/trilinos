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

#include "ReducedSpaceSQPPack/src/std/ReducedGradientStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "GeneralIterationPack/src/print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpNonsingular.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"

namespace ReducedSpaceSQPPack {

bool ReducedGradientStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	, poss_type assoc_step_poss
	)
{
	using BLAS_Cpp::trans;
	using LinAlgOpPack::V_MtV;

	rSQPAlgo    &algo   = rsqp_algo(_algo);
	rSQPState   &s      = algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	// Get iteration quantities
	IterQuantityAccess<VectorWithOpMutable>
		&Gf_iq  = s.Gf(),
		&rGf_iq = s.rGf();
	IterQuantityAccess<MatrixWithOp>
		&Z_iq = s.Z();

	// rGf = Z' * Gf
	V_MtV( &rGf_iq.set_k(0), Z_iq.get_k(0), trans, Gf_iq.get_k(0) );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||rGf||inf = "	<< rGf_iq.get_k(0).norm_inf() << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nrGf_k =\n" << rGf_iq.get_k(0);
	}

	return true;
}

void ReducedGradientStd_Step::print_step(
	const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Evaluate the reduced gradient of the objective funciton\n"
		<< L << "rGf_k = Z_k' * Gf_k\n";
}

} // end namespace ReducedSpaceSQPPack
