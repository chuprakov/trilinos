// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_Step.h

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproachTester.h"
#include "ConstrainedOptimizationPack/include/VariableBoundsTester.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Base class for evaluating a new point for the "Tailored Appraoch".
  *
  * Uses the NLPrSQPTailoredApproach interface to compute Z = [ -inv(C)*N; I ] and
  * py = -inv(C)*c explicitly.  Subclasses determine how py and Ypy are updated.
  */
class EvalNewPointTailoredApproach_Step : public EvalNewPoint_Step {
public:

	/// «std comp» members for comparision object compatible with Gc
	STANDARD_COMPOSITION_MEMBERS( NLPrSQPTailoredApproachTester, deriv_tester )

	///
	typedef ConstrainedOptimizationPack::VariableBoundsTester
		VariableBoundsTester;

	/// «std comp» Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )

	///
	enum EFDDerivTesting { FD_DEFAULT, FD_TEST, FD_NO_TEST };

	///
	/** Set how and if finite derivatives are tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDDerivTesting, fd_deriv_testing )

	/// set new_point == true by default.
	EvalNewPointTailoredApproach_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, const bounds_tester_ptr_t&	bounds_tester
		, EFDDerivTesting				fd_deriv_testing
		);

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

protected:

	///
	/** Overridden by subclass to compute py and Ypy.
	  *
	  * The matrix Y is never computed explicity.
	  *
	  * @param	D	[I]		matrix D = -inv(C)*N
	  * @param	py	[I/O]	On input py = -inv(C)*c, on output py = -inv((Gc(decomp)'*Y)*c
	  * @param	Ypy	[O]		On ouput is Y*py.
	  */
	virtual void calc_py_Ypy( const GenMatrixSlice& D, VectorSlice* py, Vector* Ypy
		, EJournalOutputLevel olevel, std::ostream& out ) = 0;

	///
	/** Overridden by subclass print how py and Ypy are computed.
	  */
	virtual void print_calc_Y_py_Ypy( std::ostream& out, const std::string& leading_str ) const = 0;

private:

	// Not defined and not to be called
	EvalNewPointTailoredApproach_Step();

};	// end class EvalNewPointTailoredApproach_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H
