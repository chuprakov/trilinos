// ////////////////////////////////////////////////////////////////////////////
// CalcLambdaIndepStd_AddedStep.h

#ifndef CALC_LAMBDA_INDEP_STD_ADDED_STEP_H
#define CALC_LAMBDA_INDEP_STD_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Calculates the lagrange multipliers for the independent constraints.
  */
class CalcLambdaIndepStd_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CalcLambdaIndepStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// CALC_LAMBDA_INDEP_STD_ADDED_STEP_H