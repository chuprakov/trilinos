// ////////////////////////////////////////////////////////////////////////////
// CalcDFromYPYZPZ_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <limits>
#include <ostream>

#include "../../include/std/CalcDFromYPYZPZ_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/print_vector_change_stats.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::CalcDFromYPYZPZ_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::V_VpV;
	using LinAlgPack::dot;
	using LinAlgPack::norm_inf;

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
	VectorWithNorms &d = s.d().set_k(0);
	V_VpV( &d.v(), s.Ypy().get_k(0)(), s.Zpz().get_k(0)() );

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		const value_type very_small = std::numeric_limits<value_type>::min();
		out << "\n(Ypy_k'*Zpz_k)/((eps+||Ypy_k||2)*(eps+||Zpz_k||2)) = "
			<< dot( s.Ypy().get_k(0)(), s.Zpz().get_k(0)() )
				/ ((very_small + s.Ypy().get_k(0).norm_2())
					* (very_small + s.Zpz().get_k(0).norm_2() ) );
		out	<< "\n||d||inf = " << d.norm_inf() << std::endl;
		ConstrainedOptimizationPack::print_vector_change_stats(
			s.x().get_k(0)(), "x", s.d().get_k(0)(), "d", out );
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nd_k = \n" << s.d().get_k(0)();
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
