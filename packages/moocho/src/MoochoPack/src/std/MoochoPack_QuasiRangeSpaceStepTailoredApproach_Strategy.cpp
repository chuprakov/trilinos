// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepTailoredApproach_Strategy.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "ReducedSpaceSQPPack/include/std/QuasiRangeSpaceStepTailoredApproach_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/rSQPAlgorithmStepNames.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproach.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"
#include "ConstrainedOptimizationPack/include/DenseIdentVertConcatMatrixSubclass.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/WorkspacePack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

bool QuasiRangeSpaceStepTailoredApproach_Strategy::solve_quasi_range_space_step(
	std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* v
  	)
{
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// Get NLP reference
#ifdef _WINDOWS
	NLPrSQPTailoredApproach
		&nlp	= dynamic_cast<NLPrSQPTailoredApproach&>(algo->nlp());
#else
	NLPrSQPTailoredApproach
		&nlp	= dyn_cast<NLPrSQPTailoredApproach>(algo->nlp());
#endif

	// Get D for Z_k = [ D; I ]
	const MatrixWithOp
		&Z_k = s->Z().get_k(0);
#ifdef _WINDOWS
	const DenseIdentVertConcatMatrixSubclass
		&cZ_k = dynamic_cast<const DenseIdentVertConcatMatrixSubclass&>(Z_k);
#else
	const DenseIdentVertConcatMatrixSubclass
		&cZ_k = dyn_cast<const DenseIdentVertConcatMatrixSubclass>(Z_k);
#endif
	const GenMatrixSlice
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
	Vector c_xo_tmp = c_xo, vy_tmp;  // This is hacked.  This sucks!
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

} // end namespace ReducedSpaceSQPPack
