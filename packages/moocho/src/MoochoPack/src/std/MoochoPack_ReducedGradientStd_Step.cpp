// ////////////////////////////////////////////////////////////////////////////
// ReducedGradientStd_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/ReducedGradientStd_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOut.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::ReducedGradientStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgOpPack::V_MtV;
	using LinAlgPack::norm_inf;

	rSQPState	&s		= rsqp_algo(_algo).rsqp_state();

	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = _algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	// rGf = Z' * Gf
	V_MtV( &s.rGf().set_k(0).v(), s.Z().get_k(0), BLAS_Cpp::trans, s.Gf().get_k(0)() );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||rGf||inf = "	<< s.rGf().get_k(0).norm_inf() << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nrGf_k =\n" << s.rGf().get_k(0)();
	}

	return true;
}

void ReducedSpaceSQPPack::ReducedGradientStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the reduced gradient of the objective funciton\n"
		<< L << "rGf_k = Z_k' * Gf_k\n";
}
