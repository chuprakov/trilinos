// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"

namespace MoochoPack {

///
/** Updates the penalty parameter for a merit function as:
 * mu_k = max( mu_km1, min_mu ).
 *
 * This class assumes the merit function iteration quantity supports
 * the interfaces <tt>MeritFuncPenaltyParam</tt> and <tt>MeritFuncNLPDirecDeriv</tt>.
 * 
 * min_mu is computed by subclasses.
 */
class MeritFunc_PenaltyParamUpdateGuts_AddedStep
	: public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

	///
	MeritFunc_PenaltyParamUpdateGuts_AddedStep(
		value_type   small_mu
		,value_type  mult_factor
		,value_type  kkt_near_sol
		);

	/** @name Overridden from AlgorithmStep */
	//@{

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss
		, IterationPack::EDoStepType type, poss_type assoc_step_poss
		, std::ostream& out, const std::string& leading_str ) const;

	//@}

	/** @name Overridden from MeritFunc_PenaltyParamUpdate_AddedStep */
	//@{

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
	
	//@}

protected:

	/** @name Pure virtual functions to be overridden by subclasses */
	//@{

	///
	/** Override to determine the mininum value of mu the penalty parameter
	 * can take on.
	 *
	 * @param  s       [in] The rSQP state object to get at useful information.
	 * @param  min_mu  [out] If min_mu(...) returns true, then this is the
	 * 					mininum value mu can take on and still have
	 * 					descent in the merit function.
	 * 					
	 * @return	Returns true if the penalty parameter should be updated
	 * 	or false if the previous mu_km1 should be used.
	 */
	virtual bool min_mu( NLPAlgoState& s, value_type* min_mu ) const = 0;

	///
	/** Override to print how min_mu calculated.
	 */
	virtual void print_min_mu_step( std::ostream& out
		, const std::string& leading_str ) const = 0;

	//@}

private:
	bool near_solution_;
	value_type small_mu_;
	value_type mult_factor_;
	value_type kkt_near_sol_;
	
};	// end class MeritFunc_PenaltyParamUpdateGuts_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
