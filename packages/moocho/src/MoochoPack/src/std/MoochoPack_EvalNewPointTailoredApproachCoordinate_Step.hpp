// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.hpp
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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"

namespace MoochoPack {

///
/** Implements "coordinate" decompostion for "Tailored Appraoch".
  *
  * Computes:<br>
  * <tt>py = py</tt><br>
  * <tt>Y = [ I; 0 ]</tt><br>
  */
class EvalNewPointTailoredApproachCoordinate_Step
	: public EvalNewPointTailoredApproach_Step
{
public:

	///
	EvalNewPointTailoredApproachCoordinate_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, const bounds_tester_ptr_t&	bounds_tester
		, EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
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
	// not defined and not to be called
	EvalNewPointTailoredApproachCoordinate_Step();

};	// end class EvalNewPointTailoredApproachCoordinate_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
