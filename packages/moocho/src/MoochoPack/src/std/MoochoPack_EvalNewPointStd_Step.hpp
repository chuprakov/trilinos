// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointStd_Step.h

#ifndef EVAL_NEW_POINT_STD_STEP_H
#define EVAL_NEW_POINT_STD_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "NLPInterfacePack/test/NLPFirstDerivativesTester.h"
#include "ConstrainedOptimizationPack/include/VariableBoundsTester.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Standard new point evaluation step class.
  *
  * This class calcualtes Gc, updates Z, and Y, and calculates Gf, c, and f in that order.
  * Also the lagrange multipliers (lambda) for the equality constraints (c) can be
  * calculated if compute_lambda(true) is called.
  */
class EvalNewPointStd_Step : public EvalNewPoint_Step {
public:

	///
	typedef NLPInterfacePack::TestingPack::NLPFirstDerivativesTester
		NLPFirstDerivativesTester;

	/// «std comp» Members for first derivative tester object
	STANDARD_COMPOSITION_MEMBERS( NLPFirstDerivativesTester, deriv_tester )

	///
	typedef ConstrainedOptimizationPack::VariableBoundsTester
		VariableBoundsTester;

	/// «std comp» Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )

	///
	/** Call to set #newx# = #new_point# which is passed to the to the first nlp calc with is Gc.
	  *
	  * This is primarily used when a new basis must be selected for variable
	  * reduction decompositions for Z, and Y.
	  *
	  * For this operation to have its effect it must be called before each call to
	  * #do_step()# and therefore the effect does not presist between calls of
	  * #do_step()#.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, new_point )

	///
	enum EFDDerivTesting { FD_DEFAULT, FD_TEST, FD_NO_TEST };

	///
	/** Set how and if finite derivatives are tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDDerivTesting, fd_deriv_testing )

	/// set new_point == true by default.
	EvalNewPointStd_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, const bounds_tester_ptr_t&	bounds_tester
		, EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
		);

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:

	// Not defined and not to be called
	EvalNewPointStd_Step();

};	// end class EvalNewPointStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_STD_STEP_H
