// ////////////////////////////////////////////////////////////////////////////
// NumFixedDepIndep_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/NumFixedDepIndep_AddedStep.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"

bool ReducedSpaceSQPPack::NumFixedDepIndep_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	Range1D		indep	= s.con_indep();

	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	if( s.nu().updated_k(0) && s.nu().get_k(0).nz() ) {
		const Range1D
			dep		= s.var_dep(),
			indep	= s.var_indep();
		const SpVector &nu_k	= s.nu().get_k(0);
		size_type fixed_dep = 0, fixed_indep = 0;
		for( SpVector::const_iterator itr = nu_k.begin(); itr != nu_k.end(); ++itr ) {
			if( dep.in_range( itr->indice() ) )
				fixed_dep++;
			else if( indep.in_range( itr->indice() ) )
				fixed_indep++;
			else
				assert(0);	// should never happen
		}
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nnum_dep_fixed = "		<< fixed_dep
				<< "\nnum_indep_fixed = "	<< fixed_indep << std::endl;
		}
	}
	else {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nnu not calculated for the kth iteration\n";
		}
	}

	return true;
}

void ReducedSpaceSQPPack::NumFixedDepIndep_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Counts the number of fixed variables from "
				"the dependent and independent sets\n";
}
