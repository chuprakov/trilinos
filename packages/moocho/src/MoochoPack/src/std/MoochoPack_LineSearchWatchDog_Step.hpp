// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_LineSearchWatchDog_Step.hpp
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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_H
#define LINE_SEARCH_WATCH_DOG_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Implements watchdog line search.
  *
  * The watchdog proceedure will only be considered when :
  * opt_kkt_err_k <= opt_kkt_err_threshold && feas_kkt_err_k <= feas_kkt_err_threshold
  * (see step listing).  The default behavior is to never use the watchdog procedure.
  */
class LineSearchWatchDog_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for direct_line_search
	STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search)

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	/** <<std member comp>> members for the armijo fractional reduction parameter.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta )

	///
	/** <<std member comp>> members for the threshold for opt_kkt_err before
	  * the watchdog procedure should kick-in.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_kkt_err_threshold )
	
	///
	/** <<std member comp>> members for the threshold for feas_kkt_err before
	  * the watchdog procedure should kick-in.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_kkt_err_threshold )

	///
	LineSearchWatchDog_Step(
			  const direct_line_search_ptr_t&	direct_line_search		= 0
			, const merit_func_ptr_t&			merit_func				= 0
			, value_type						eta						= 1e-4
			, value_type						opt_kkt_err_threshold 	= 1e-1
			, value_type						feas_kkt_err_threshold 	= 1e-3
			);

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	int					watch_k_;	// >= 0 means that we are using the watchdog.
	DVector				xo_;
	value_type			fo_;
	value_type			nrm_co_;
	DVector				do_;
	value_type			phio_;
	value_type			Dphio_;
	value_type			phiop1_;
	
};	// end class LineSearchWatchDog_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_WATCH_DOG_STEP_H
