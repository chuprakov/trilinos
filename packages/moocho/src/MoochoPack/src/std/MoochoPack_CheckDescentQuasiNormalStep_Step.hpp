// ////////////////////////////////////////////////////////////////////
// CheckDescentQuasiNormalStep_Step.hpp
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

#ifndef CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H
#define CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H

#include "MoochoPack/src/MoochoPackTypes.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "NLPInterfacePack/src/abstract/tools/CalcFiniteDiffProd.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Checks for descent in the decomposed equality constraints with respect to the
 * range space step <tt>Ypy</tt> using finite differences.
 *
 * This step class checks for descent in the feasibility measure <tt>q(x) = 1/2 * cd(x)'*cd(x) <: REAL</tt>
 * of the decomposed equality constraints <tt>cd(x) = c(equ_decomp)(x)</tt> with respect to the range space
 * step <tt>Ypy_k</tt>.  The gradient of this feasibility measure is:
 \verbatim

 grad(q(x),x) = grad(cd(x),x) * cd(x)
 \endverbatim
 * Therefore, we can determine if we have descent by checking <tt>grad(q(x),x)'*Ypy_k = cd(x)'*grad(cd(x),x)'*Ypy_k< 0</tt>.
 * The product <tt>grad(c(x),x)'*Ypy_k</tt> is approximated with finite differences using the class
 * <tt>MoochoPack::CalcFiniteDiffProd</tt>.
 */
class CheckDescentQuasiNormalStep_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/// Set the object that will compute the finite difference products.
	STANDARD_COMPOSITION_MEMBERS( CalcFiniteDiffProd, calc_fd_prod )

	///
	/** Constructor
	 */
	CheckDescentQuasiNormalStep_Step(
		const calc_fd_prod_ptr_t&   calc_fd_prod
		);

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

};	// end class CheckDescentQuasiNormalStep_Step

}	// end namespace MoochoPack 

#endif	// CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H
