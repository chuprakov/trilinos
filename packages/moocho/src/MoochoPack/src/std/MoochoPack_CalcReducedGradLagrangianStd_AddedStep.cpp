// ////////////////////////////////////////////////////////////////////////////
// CalcReducedGradLagrangianStd_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include "../../include/std/CalcReducedGradLagrangianStd_AddedStep.h"
#include "../../include/rSQPAlgoContainer.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorOut.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::CalcReducedGradLagrangianStd_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using BLAS_Cpp::trans;
	using LinAlgPack::norm_inf;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::Vp_MtV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Calculate rGL = Z' * Gf + Z' * nu + V' * lambda(dep)

	VectorWithNorms &rGL = s.rGL().set_k(0);

	if( s.nu().updated_k(0) ) {
		// Compute rGL = Z'*(Gf + nu) to reduce the effect of roundoff in this
		// catastropic cancelation.
		Vector tmp;	// tmp = Gf + nu
		tmp = s.Gf().get_k(0)();
		Vp_V( &tmp(), s.nu().get_k(0)() );
		V_MtV(	&rGL.v(), s.Z().get_k(0), trans, tmp() );
	}
	else {
		rGL.v() = s.rGf().get_k(0)();
	}

//	int stupid = 1;	// Put is to avoid an internal compiler error in Release mode.
							
	// rGL += V' * lambda(dep)					
	if( algo.nlp().r() < algo.nlp().m() )
		Vp_MtV( &rGL.v()(), s.V().get_k(0), trans, s.lambda().get_k(0).v()(s.con_undecomp()) );	

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||rGL||inf = " << rGL.norm_inf() << "\n";
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nrGL_k = \n" << rGL.v()();
	}

	return true;
}

void ReducedSpaceSQPPack::CalcReducedGradLagrangianStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the reduced gradient of the Lagrangian\n"
		<< L << "if nu_k is updated then\n"
		<< L << "    rGL_k = Z_k' * (Gf_k + nu_k) + V_k' * lambda_k(undecomp)\n"
		<< L << "else\n"
		<< L << "    rGL_k = rGf_k + V_k' * lambda_k(undecomp)\n"
		<< L << "end\n";
}
