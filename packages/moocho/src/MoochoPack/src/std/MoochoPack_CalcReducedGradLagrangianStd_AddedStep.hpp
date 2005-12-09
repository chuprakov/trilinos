// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_CalcReducedGradLagrangianStd_AddedStep.hpp
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

#ifndef CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
#define CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

///
/** Calculates the reduced gradient of the Lagrangian
 * <tt>rGL = rGf + Z' * nu + GcUP' * lambda(equ_undecomp) + GhUP' * lambdaI(inequ_undecomp)</tt> 
 */
class CalcReducedGradLagrangianStd_AddedStep
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CalcReducedGradLagrangianStd_AddedStep

}	// end namespace MoochoPack 

#endif	// CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
