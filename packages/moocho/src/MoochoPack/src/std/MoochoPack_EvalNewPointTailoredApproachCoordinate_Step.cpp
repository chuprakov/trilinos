// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachCoordinate_Step.cpp

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachCoordinate_Step.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace ReducedSpaceSQPPack {

EvalNewPointTailoredApproachCoordinate_Step::EvalNewPointTailoredApproachCoordinate_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, EFDDerivTesting				fd_deriv_testing
		)
	:
		EvalNewPointTailoredApproach_Step(deriv_tester,fd_deriv_testing)
{}
// protected

void EvalNewPointTailoredApproachCoordinate_Step::calc_py_Ypy(
	  const GenMatrixSlice& D, VectorSlice* py, Vector* Ypy
	, EJournalOutputLevel olevel, std::ostream& out
	)
{
	const size_type
		n = D.rows()+D.cols(),
		r = D.cols();
	// py is not altered
	Ypy->resize(n);
	(*Ypy)(1,r) = *py;
	(*Ypy)(r+1,n) = 0.0;
}

void EvalNewPointTailoredApproachCoordinate_Step::print_calc_Y_py_Ypy(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Coordinate decomposition\n"
		<< L << "py_k = py_k\n"
		<< L << "Y = [ I ; 0 ] <: R^(n x m)   [Not computed explicity]\n"
		<< L << "Ypy_k = Y * py_k\n"
		;
}

}	// end namespace ReducedSpaceSQPPack 
