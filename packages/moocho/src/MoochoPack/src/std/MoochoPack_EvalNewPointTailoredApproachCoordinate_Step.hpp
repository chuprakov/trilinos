// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachCoordinate_Step.h

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Implements "coordinate" decompostion for "Tailored Appraoch".
  *
  * Computes:\\
  * py = py \\
  * Ypy = [ py; 0 ] \\
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

	// ///////////////////////////////
	// Overridden

	///
	void calc_py_Ypy( const GenMatrixSlice& D, VectorSlice* py, Vector* Ypy
		, EJournalOutputLevel olevel, std::ostream& out );

	///
	void print_calc_Y_py_Ypy( std::ostream& out, const std::string& leading_str ) const;

private:
	// not defined and not to be called
	EvalNewPointTailoredApproachCoordinate_Step();

};	// end class EvalNewPointTailoredApproachCoordinate_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
