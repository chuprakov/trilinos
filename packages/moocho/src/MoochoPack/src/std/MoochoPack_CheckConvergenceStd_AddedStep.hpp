// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.h
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

#ifndef CHECK_CONVERGENCE_STD_ADDEDSTEP_H
#define CHECK_CONVERGENCE_STD_ADDEDSTEP_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/AlgorithmStep.h"
#include "StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Check for convergence.
  */
class CheckConvergenceStd_AddedStep
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
{
public:

	///
	enum EOptErrorCheck { OPT_ERROR_REDUCED_GRADIENT_LAGR, OPT_ERROR_GRADIENT_LAGR };

	///
	/** <<std member comp>> members for whether to check the reduced
	  * or full gradient of the Lagrangian. 
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EOptErrorCheck, opt_error_check )

	///
	enum EScaleKKTErrorBy { SCALE_BY_ONE, SCALE_BY_NORM_2_X, SCALE_BY_NORM_INF_X };

	///
	/** <<std member comp>> members for whether the optimality conditions
	  * should be scaled by the gradient of the objective or not.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EScaleKKTErrorBy, scale_kkt_error_by )

	///
	/** <<std member comp>> members for whether the optimality conditions
	  * should be scaled by the 
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, scale_opt_error_by_Gf )

	///
	CheckConvergenceStd_AddedStep(
		  EOptErrorCheck opt_error_check		= OPT_ERROR_REDUCED_GRADIENT_LAGR
		, EScaleKKTErrorBy scale_kkt_error_by	= SCALE_BY_ONE
		, bool scale_opt_error_by_Gf 			= true
		);

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

};	// end class CheckConvergenceStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// CHECK_CONVERGENCE_STD_ADDEDSTEP_H
