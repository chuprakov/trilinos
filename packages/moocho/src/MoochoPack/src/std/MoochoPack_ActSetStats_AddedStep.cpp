// ////////////////////////////////////////////////////////////////////////////
// ActSetStats_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/ActSetStats_AddedStep.h"
#include "../../include/std/active_set_change.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"

bool ReducedSpaceSQPPack::ActSetStats_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( s.nu().updated_k(0) ) {
		size_type num_active = 0, num_adds = 0, num_drops = 0;
		const SpVector &nu_k	= s.nu().get_k(0);
		num_active = nu_k.nz();
		if( s.nu().updated_k(-1) ) {
			const SpVector &nu_km1	= s.nu().get_k(-1);
			active_set_change( nu_k(), nu_km1(), olevel, &num_adds, &num_drops, &out );
			act_set_stats_(s).set_k(0).set_stats(num_active,num_adds,num_drops);
		}
		else {
			act_set_stats_(s).set_k(0).set_stats( num_active, ActSetStats::NOT_KNOWN
				,ActSetStats::NOT_KNOWN );
		}
	}
	else {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nnu not calculated for the kth iteration\n";
		}
	}

	return true;
}

void ReducedSpaceSQPPack::ActSetStats_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Updates active set statistics for changes from the last iteration\n"
		<< L << "Given nu_km1 and nu_k update act_set_stats_k (num_active, num_adds, num_drops)\n";
}
