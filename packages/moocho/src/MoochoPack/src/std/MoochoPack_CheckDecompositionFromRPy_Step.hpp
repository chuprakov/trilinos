// ////////////////////////////////////////////////////////////////////////////
// CheckDecompositionFromRPy_Step.hpp
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

#ifndef CHECK_DECOMPOSITION_FROM_R_PY_STEP_H
#define CHECK_DECOMPOSITION_FROM_R_PY_STEP_H

#include "NewDecompositionSelection_Strategy.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Check if the decomposition is going singular and if it is select a new decomposition.
 *
 * This steps checks if the decomposition is going singular if the computation
 * for the range space step looks like it is becomming more inaccurate.
 *
 * In particular we want to know how cond(C) is changing.  To do this we
 * will monitor increases in the error in solving the equation:
 \verbatim

   R*py + c(equ_decomp) = 0
 \endverbatim
 * Therefore we will check for increases in the ratio:
 \verbatim

   ratio = ||R*py + c(equ_decomp)||inf / ||c(equ_decomp)||inf
 \endverbatim
 *
 * If this ratio goes up dramatically, then this is a tail tell
 * sign that <tt>R</tt> is becomming illconditioned.  See the algorithm
 * printout for a more detailed description what what is going on.
 */
class CheckDecompositionFromRPy_Step
	: public IterationPack::AlgorithmStep // doxygen needs full name
{
public:

	/// <<std comp>> members for Decomposition Select Strategy object.
	STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy )

	///
	/** Set the maximum change in the relative error in the range space
     * step before the selection of a new decomposition is triggered.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_decomposition_cond_change_frac )

	///
	CheckDecompositionFromRPy_Step(
		 const new_decomp_strategy_ptr_t    &new_decomp_strategy
		,value_type                         max_decomposition_cond_change_frac = 1e+4
		);

	/// Call the reset initialization of all defaults.
	void reset();

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:

	value_type	beta_min_;

	// Not defined and not to be called
	CheckDecompositionFromRPy_Step();

};	// end class CheckDecompositionFromRPy_Step

}	// end namespace MoochoPack 

#endif	// CHECK_DECOMPOSITION_FROM_R_PY_STEP_H
