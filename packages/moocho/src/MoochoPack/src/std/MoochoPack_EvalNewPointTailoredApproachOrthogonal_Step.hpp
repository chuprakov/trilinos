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

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"
#include "LinAlgPack/include/GenMatrixClass.h"

namespace ReducedSpaceSQPPack {

///
/** Implements "orthogonal" decompostion for "Tailored Appraoch".
  *
  * Computes:\\
  * py = inv(I + D*D') * py \\
  * Ypy = [ py; -D'*py ] \\
  */
class EvalNewPointTailoredApproachOrthogonal_Step
	: public EvalNewPointTailoredApproach_Step
{
public:

	///
	EvalNewPointTailoredApproachOrthogonal_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, const bounds_tester_ptr_t&	bounds_tester
		, EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
		);

protected:

	// ///////////////////////////////
	// Overridden

	///
	void calc_py_Ypy( const GenMatrixSlice& D, VectorSlice* py, VectorSlice* Ypy
		, EJournalOutputLevel olevel, std::ostream& out );
	///
	void recalc_py_Ypy( const GenMatrixSlice& D, VectorSlice* py, VectorSlice* Ypy
		, EJournalOutputLevel olevel, std::ostream& out );
	///
	void print_calc_Y_py_Ypy( std::ostream& out, const std::string& leading_str ) const;

private:

	// ///////////////////////////////
	// Private data members

	GenMatrix   SL_;

	// //////////////////////////////
	// Private member functions

	// not defined and not to be called
	EvalNewPointTailoredApproachOrthogonal_Step();

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
