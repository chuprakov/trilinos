// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachOrthogonal_Step.hpp
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

#include "EvalNewPointTailoredApproach_Step.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

///
/** Implements "orthogonal" decompostion for "Tailored Appraoch".
 *
 * Computes:
 \verbatim
 py = inv(I + D*D') * py
 Y  = [ I; -D' ]
 Uy = ???
 \endverbatim
 */
class EvalNewPointTailoredApproachOrthogonal_Step
	: public EvalNewPointTailoredApproach_Step
{
public:

	///
	EvalNewPointTailoredApproachOrthogonal_Step(
		const deriv_tester_ptr_t                &deriv_tester
		,const bounds_tester_ptr_t              &bounds_tester
		,EFDDerivTesting                        fd_deriv_testing = FD_DEFAULT
		);

protected:

	/** @name Overridden from EvalNewPointTailoredApproach_Step */
	//@{

	///
	void uninitialize_Y_Uy(
		MatrixOp         *Y
		,MatrixOp        *Uy
		);
	///
	void calc_py_Y_Uy(
		const NLPDirect       &nlp
		,const D_ptr_t        &D
		,VectorMutable        *py
		,MatrixOp             *Y
		,MatrixOp             *Uy
		,EJournalOutputLevel  olevel
		,std::ostream         &out
		);
	///
	void recalc_py(
		const MatrixOp           &D
		,VectorMutable           *py
		,EJournalOutputLevel     olevel
		,std::ostream            &out
		);
	///
	void print_calc_py_Y_Uy(
		std::ostream& out, const std::string& leading_str
		) const;

	//@}

private:

	// ///////////////////////////////
	// Private types

	///
	typedef Teuchos::RefCountPtr<MatrixSymOpNonsing>  S_ptr_t;

	// ///////////////////////////////
	// Private data members

	S_ptr_t   S_ptr_;

	// //////////////////////////////
	// Private member functions

	// not defined and not to be called
	EvalNewPointTailoredApproachOrthogonal_Step();

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
