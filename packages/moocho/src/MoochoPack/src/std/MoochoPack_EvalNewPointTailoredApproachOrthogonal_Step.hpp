// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachOrthogonal_Step.h
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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H

#include "EvalNewPointTailoredApproach_Step.h"
#include "ConstrainedOptimizationPack/include/VarReductOrthog_Strategy.h"
#include "StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Implements "orthogonal" decompostion for "Tailored Appraoch".
 *
 * Computes:
 \verbatim
 py = inv(I + D*D') * py
 Y  = [ I; -D' ]
 Uy = ???
 Vy = ???
 \endverbatim
 */
class EvalNewPointTailoredApproachOrthogonal_Step
	: public EvalNewPointTailoredApproach_Step
{
public:

	///
	STANDARD_CONST_COMPOSITION_MEMBERS( VarReductOrthog_Strategy, var_reduct_orthog_strategy )

	///
	EvalNewPointTailoredApproachOrthogonal_Step(
		const var_reduct_orthog_strategy_ptr_t  &var_reduct_orthog_strategy
		,const deriv_tester_ptr_t               &deriv_tester
		,const bounds_tester_ptr_t              &bounds_tester
		,EFDDerivTesting                        fd_deriv_testing = FD_DEFAULT
		);

protected:

	/** @name Overridden from EvalNewPointTailoredApproach_Step */
	//@{

	///
	void uninitialize_Y_Uv_Uy(
		MatrixWithOp         *Y
		,MatrixWithOp        *Uy
		,MatrixWithOp        *Vy
		);
	///
	void calc_py_Y_Uy_Vy(
		const NLPFirstOrderDirect   &nlp
		,const D_ptr_t              &D
		,VectorWithOpMutable        *py
		,MatrixWithOp               *Y
		,MatrixWithOp               *Uy
		,MatrixWithOp               *Vy
		,EJournalOutputLevel        olevel
		,std::ostream               &out
		);
	///
	void recalc_py(
		const MatrixWithOp       &D
		,VectorWithOpMutable     *py
		,EJournalOutputLevel     olevel
		,std::ostream            &out
		);
	///
	void print_calc_py_Y_Uy_Vy(
		std::ostream& out, const std::string& leading_str
		) const;

	//@}

private:

	// ///////////////////////////////
	// Private types

	///
	typedef ReferenceCountingPack::ref_count_ptr<MatrixSymWithOpNonsingular>  S_ptr_t;

	// ///////////////////////////////
	// Private data members

	S_ptr_t   S_ptr_;

	// //////////////////////////////
	// Private member functions

	// not defined and not to be called
	EvalNewPointTailoredApproachOrthogonal_Step();

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
