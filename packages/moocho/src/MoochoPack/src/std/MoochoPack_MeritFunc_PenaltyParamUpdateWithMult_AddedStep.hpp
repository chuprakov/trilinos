// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateWithMult_AddedStep.h

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULT_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULT_ADDED_STEP_H

#include "MeritFunc_PenaltyParamUpdate_AddedStep.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Updates the penalty parameter for a merit function as:
  * mu_k = max( mu_km1, ||lambda||inf ).
  *
  * This class assumes the merit function supports the interfaces
  * MeritFuncPenaltyParam and MeritFuncNLPDirecDeriv.
  */
class MeritFunc_PenaltyParamUpdateWithMult_AddedStep
	: public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	MeritFunc_PenaltyParamUpdateWithMult_AddedStep(const merit_func_ptr_t& merit_func
		, value_type small_mu = 1e-6, value_type mult_factor = 1e-4
		, value_type kkt_near_sol = 1.0 );

	// ///////////////////////////////
	// Overridden from AlgorithmStep

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss
		, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
		, std::ostream& out, const std::string& leading_str ) const;

	// //////////////////////////////////////////////////////////////////////
	// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

	///
	void small_mu( value_type small_mu );

	///
	value_type small_mu() const;

	///
	void mult_factor( value_type mult_factor );

	///
	value_type mult_factor() const;

	///
	void kkt_near_sol( value_type kkt_near_sol );

	///
	value_type kkt_near_sol() const;

private:
	bool near_solution_;
	value_type small_mu_;
	value_type mult_factor_;
	value_type kkt_near_sol_;
	
};	// end class MeritFunc_PenaltyParamUpdateWithMult_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULT_ADDED_STEP_H