// ////////////////////////////////////////////////////////////////////////////
// CheckSkipBFGSUpdateStd_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/CheckSkipBFGSUpdateStd_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "../../include/std/quasi_newton_stats.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

bool ReducedSpaceSQPPack::CheckSkipBFGSUpdateStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_2;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( s.Ypy().updated_k(-1) && s.Zpz().updated_k(-1)
			&& s.rGL().updated_k(-1) && s.rHL().updated_k(-1)
			&& s.c().updated_k(-1) )
	{
		// The information exists for the update so determine
		// if we are in the region to perform the BFGS update.
		
		// Check if we are to skip the update for this iteration
		
		const value_type
			nrm_rGL_km1 = s.rGL().get_k(-1).norm_2(),
			nrm_c_km1	= s.c().get_k(-1).norm_2(),
			nrm_Zpz_km1	= s.Zpz().get_k(-1).norm_2(),
			nrm_Ypy_km1	= s.Ypy().get_k(-1).norm_2();

		// ratio = (10.0 / sqrt(||rGL_km1|| + ||c_km1||)) * ( ||Zpz_km1|| / ||Ypy_km1|| )
		value_type
			ratio = ( 10.0 / ::sqrt( nrm_rGL_km1 + nrm_c_km1 ) )
				* ( nrm_Zpz_km1 / nrm_Ypy_km1 );

		// If ratio > 1.0 then skip the update
		const bool skip_update = ratio < 1.0;

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out
				<< "ratio = (10/sqrt(||rGL_km1||+||c_km1||))*(||Zpz_km1||/||Ypy_km1||)\n"
				<< "      = (10/sqrt("<<nrm_rGL_km1<<"+"<<nrm_c_km1<<"))\n"
				<< "        * ("<<nrm_Zpz_km1<<"/"<<nrm_Ypy_km1<<")\n"
				<< "      = " << ratio << std::endl
				<< "ratio " << (skip_update ? '<' : '>' ) << " 1\n"
				<< (skip_update
						? "Skipping BFGS update ...\n"
						: "Perform BFGS update if you can ...\n"  );
		}	

		if(	skip_update ) {

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out
					<< "\nrHL_k = rHL_km1\n";
			}

			const MatrixWithOp &rHL_km1 = s.rHL().get_k(-1);
			s.rHL().set_k(0) = rHL_km1;
			quasi_newton_stats(s).set_k(0).set_updated_stats(
				QuasiNewtonStats::SKIPED );

			if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
				s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
			}
		}

	}

	return true;
}

void ReducedSpaceSQPPack::CheckSkipBFGSUpdateStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Check if we should do the BFGS update\n"
		<< L << "if Ypy, Zpz, rGL, rHL and c are updated for the k-1 iteration\n"
		<< L << "    *** Check if we are in the proper region\n"
		<< L << "    ratio = ( 10 / sqrt( norm(rGL_km1,2) + norm(c_km1,2) ) )\n"
		<< L << "             * ( norm(Zpz_km1,2) / norm(Ypy_km1,2) )\n"
		<< L << "    if ratio < 1 then \n"
		<< L << "        rHL_k = rHL_km1\n"
		<< L << "    end\n"
		<< L << "end\n";
}
