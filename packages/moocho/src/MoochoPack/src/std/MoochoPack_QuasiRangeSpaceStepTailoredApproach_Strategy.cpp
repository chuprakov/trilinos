// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepTailoredApproach_Strategy.cpp
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

#include "MoochoPack_QuasiRangeSpaceStepTailoredApproach_Strategy.hpp"
#include "MoochoPack_MoochoAlgorithmStepNames.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "MoochoPack/src/NLPrSQPTailoredApproach.h"
#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "ConstrainedOptPack/src/DenseIdentVertConcatMatrixSubclass.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "MiWorkspacePack.h"
#include "Midynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

namespace MoochoPack {

bool QuasiRangeSpaceStepTailoredApproach_Strategy::solve_quasi_range_space_step(
	std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
	,const DVectorSlice& xo, const DVectorSlice& c_xo, DVectorSlice* v
  	)
{
	using Teuchos::dyn_cast;
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	// Get NLP reference
#ifdef _WINDOWS
	NLPrSQPTailoredApproach
		&nlp	= dynamic_cast<NLPrSQPTailoredApproach&>(algo->nlp());
#else
	NLPrSQPTailoredApproach
		&nlp	= dyn_cast<NLPrSQPTailoredApproach>(algo->nlp());
#endif

	// Get D for Z_k = [ D; I ]
	const MatrixOp
		&Z_k = s->Z().get_k(0);
#ifdef _WINDOWS
	const DenseIdentVertConcatMatrixSubclass
		&cZ_k = dynamic_cast<const DenseIdentVertConcatMatrixSubclass&>(Z_k);
#else
	const DenseIdentVertConcatMatrixSubclass
		&cZ_k = dyn_cast<const DenseIdentVertConcatMatrixSubclass>(Z_k);
#endif
	const DMatrixSlice
		D = cZ_k.m().D();

	// Get reference to EvalNewPoint step
#ifdef _WINDOWS
	EvalNewPointTailoredApproach_Step
		&eval_tailored = dynamic_cast<EvalNewPointTailoredApproach_Step&>(
			*algo->get_step(algo->get_step_poss(EvalNewPoint_name)));
#else
	EvalNewPointTailoredApproach_Step
		&eval_tailored = dyn_cast<EvalNewPointTailoredApproach_Step>(
			*algo->get_step(algo->get_step_poss(EvalNewPoint_name)));
#endif

	// Compute an approximate newton step for constriants wy
	DVector c_xo_tmp = c_xo, vy_tmp;  // This is hacked.  This sucks!
	nlp.calc_semi_newton_step(xo,&c_xo_tmp,false,&vy_tmp);
		
	// Compute wy, Ywy
	eval_tailored.recalc_py_Ypy(D,&vy_tmp(),v,olevel,out);

	return true;
}

void QuasiRangeSpaceStepTailoredApproach_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out << L << "*** Compute the approximate range space step by calling on the \"Tailored Approach\" NLP interface:\n"
		<< L << "Compute vy s.t. ||Gc_k'*Y_k*vy + c_xo|| << ||c_xo|| (nlp.calc_semi_newton_step(...))\n"
		<< L << "update vy and compute v = Yvy from EvalNewPointTailoredApproach_Step::recalc_py_Ypy(...)\n";
		;
}

} // end namespace MoochoPack
