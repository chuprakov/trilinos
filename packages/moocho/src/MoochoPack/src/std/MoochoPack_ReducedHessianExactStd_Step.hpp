// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianExactStd_Step.h

#ifndef REDUCED_HESSIAN_EXACT_STD_STEP_H
#define REDUCED_HESSIAN_EXACT_STD_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"
#include "../rSQPAlgo.h"

namespace ReducedSpaceSQPPack {

///
/** Computes the exact reduced Hessian rHL_k = Z_k' * HL_k * Z_k
  */
class ReducedHessianExactStd_Step : public ReducedHessian_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class ReducedGradientStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// REDUCED_HESSIAN_EXACT_STD_STEP_H
