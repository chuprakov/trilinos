// ////////////////////////////////////////////////////////////////////////////
// NullSpaceStepWithoutBounds_Step.h
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

#ifndef NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
#define NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/AlgorithmStep.h"

namespace ReducedSpaceSQPPack {

///
/** Solves the unconstrained QP subproblem: <tt>min  qp_grad' * pz + (1/2) * pz' * rHL * pz</tt>.
  *
  * The solution to this system is just:<br>
  * <tt>pz = inv(rHL) *qp_grad</tt>.
  *
  * If use_qp_correc is false then:<br>
  *   <tt>qp_grad = rGf</tt>
  * else<br>
  *   <tt>qp_grad = rGf + zeta * ZtHLYpy.<br>
  *
  * Then <tt>Zpz = Z * pz</tt>
  */
class NullSpaceStepWithoutBounds_Step
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

};	// end class NullSpaceStepWithoutBounds_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
