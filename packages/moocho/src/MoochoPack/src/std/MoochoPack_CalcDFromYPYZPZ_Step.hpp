// ////////////////////////////////////////////////////////////////////////////
// CalcDFromYPYZPZ_Step.h

#ifndef CALC_D_FROM_YPY_ZPZ_STEP_H
#define CALC_D_FROM_YPY_ZPZ_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"

namespace ReducedSpaceSQPPack {

///
/** Calculates d = Ypy + Zpz
  */
class CalcDFromYPYZPZ_Step : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CalcDFromYPYZPZ_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// CALC_D_FROM_YPY_ZPZ_STEP_H