// ////////////////////////////////////////////////////////////////////////////
// NumFixedDepIndep_AddedStep.h

#ifndef NUM_FIXED_DEP_INDEP_ADDED_STEP_H
#define NUM_FIXED_DEP_INDEP_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Computes and outputs the number of fixed variables from the dependent
  * and independent set..
  *
  * These statistics only make sense for variable reduction decompositions.
  */
class NumFixedDepIndep_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class NumFixedDepIndep_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// NUM_FIXED_DEP_INDEP_ADDED_STEP_H