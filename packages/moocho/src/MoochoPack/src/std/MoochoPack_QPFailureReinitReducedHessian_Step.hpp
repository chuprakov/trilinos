// ////////////////////////////////////////////////////////////////////////////
// QPFailureReinitReducedHessian_Step.h

#ifndef QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H
#define QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Directs the algorithm to reinitalize the reduced Hessian on the event
 * of a QP failure.
 *
 * If the delegated Step object throws a QPFailure exception
 * then this Step object wipes out all reduced Hessian info rHL
 * for the current and previous iterations and then directs the algorithm
 * back to the ReducedHessian step (see the printed step description).
 */
class QPFailureReinitReducedHessian_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for LineSearch object.
	STANDARD_COMPOSITION_MEMBERS( rSQPAlgo_Step, null_space_step )

	///
	QPFailureReinitReducedHessian_Step( const null_space_step_ptr_t& null_space_step );

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	int last_qp_failure_k_;

	// not defined and not to be called
	QPFailureReinitReducedHessian_Step();
	QPFailureReinitReducedHessian_Step(const QPFailureReinitReducedHessian_Step&);
	QPFailureReinitReducedHessian_Step& operator=(const QPFailureReinitReducedHessian_Step&);
	
};	// end class QPFailureReinitReducedHessian_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H
