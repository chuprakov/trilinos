// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H

#include "MeritFunc_PenaltyParamUpdateGuts_AddedStep.h"

namespace ReducedSpaceSQPPack {

///
/** Specializes the update of the penalty parameter for a merit function as:
  * min_mu = |(Gf_k+nu_k)'* Ypy_k| / ||c_k||1.
  */
class MeritFunc_PenaltyParamUpdateMultFree_AddedStep
	: public MeritFunc_PenaltyParamUpdateGuts_AddedStep
{
public:

	///
	MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
		const merit_func_ptr_t& merit_func
		, value_type small_mu = 1e-6
		, value_type mult_factor = 1e-4
		, value_type kkt_near_sol = 1.0
		);

protected:

	// /////////////////////////////////////////////////////////////
	// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

	///
	bool min_mu( rSQPState& s, value_type* min_mu ) const;

	///
	void print_min_mu_step( std::ostream& out
		, const std::string& leading_str ) const;
	
};	// end class MeritFunc_PenaltyParamUpdateMultFree_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
