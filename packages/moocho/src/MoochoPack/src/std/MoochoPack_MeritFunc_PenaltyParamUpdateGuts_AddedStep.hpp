// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateGuts_AddedStep.h

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H

#include "MeritFunc_PenaltyParamUpdate_AddedStep.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Updates the penalty parameter for a merit function as:
  * mu_k = max( mu_km1, min_mu ).
  *
  * This class assumes the merit function supports the interfaces
  * MeritFuncPenaltyParam and MeritFuncNLPDirecDeriv.
  * 
  * min_mu is computed by subclasses.
  */
class MeritFunc_PenaltyParamUpdateGuts_AddedStep
	: public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	MeritFunc_PenaltyParamUpdateGuts_AddedStep(
		const merit_func_ptr_t& merit_func
		, value_type small_mu
		, value_type mult_factor
		, value_type kkt_near_sol
		);

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

protected:

	// ///////////////////////////////////////////////////
	// Template methods to be overridden by subclasses

	///
	/** Override to determine the mininum value of mu the penalty parameter
	  * can take on.
	  *
	  *	@param	s		[I]	The rSQP state object to get at useful
	  *						information.
	  * @param	min_mu	[O]	If min_mu(...) returns true, then this is the
	  * 					mininum value mu can take on and still have
	  * 					descent in the merit function.
	  * 					
	  * @return	Returns true if the penalty parameter should be updated
	  * 	or false if the previous mu_km1 should be used.
	  */
	virtual bool min_mu( rSQPState& s, value_type* min_mu ) const = 0;

	///
	/** Override to print how min_mu calculated.
	  */
	virtual void print_min_mu_step( std::ostream& out
		, const std::string& leading_str ) const = 0;

private:
	bool near_solution_;
	value_type small_mu_;
	value_type mult_factor_;
	value_type kkt_near_sol_;
	
};	// end class MeritFunc_PenaltyParamUpdateGuts_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
