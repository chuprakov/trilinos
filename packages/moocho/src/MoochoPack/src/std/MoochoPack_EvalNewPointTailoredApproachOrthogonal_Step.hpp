// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachOrthogonal_Step.h

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"

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
		  const deriv_tester_ptr_t& 	deriv_tester = 0
		, EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
		);

protected:

	// ///////////////////////////////
	// Overridden

	///
	void calc_py_Ypy( const GenMatrixSlice& D, VectorSlice* py, Vector* Ypy
		, EJournalOutputLevel olevel, std::ostream& out );

	///
	void print_calc_Y_py_Ypy( std::ostream& out, const std::string& leading_str ) const;

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
