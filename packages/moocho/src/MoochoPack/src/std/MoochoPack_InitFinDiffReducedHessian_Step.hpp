// ////////////////////////////////////////////////////////////////////////////
// InitFinDiffReducedHessian_Step.h

#ifndef INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H
#define INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H

#include "ReducedSpaceSQPPack/include/rSQPAlgo_StepBaseClasses.h"
#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Initializes the reduced hessian using a single finite difference
  * along the null space of the constraints.
  *
  * A single finite difference correction is computed along:\\
  *
  * x_fd = x_k + u * Z * e
  *
  * The step length is set to u = step_scale / ||Z*e||inf.  The
  * step length is cut back if the point x_fd is outside the
  * relaxed variable bounds.
  *
  * The finite difference is then computed as:
  *
  * rGf_fd = ( Z_k' * g(x_k + u * Z*e) - rGf_k ) / u
  *
  * The diagonal terms of the reduced hessian are then set
  * as:
  *
  * diag(i) = max( ||rGf_fd||inf , smallest_ele )  if initialization_method == SCALE_IDENTITY\\
  * diag(i) = max( rGf_fd(i)     , smallest_ele )  if initialization_method == SCALE_DIAGONAL\\
  * diag(i) = max( abs(rGf_fd(i)), smallest_ele )  if initialization_method == SCALE_DIAGONAL_ABS\\
  *
  * Where:
  *
  * smallest_ele = max( ||rGf_fd||inf / max_cond , min_diag )
  *
  * Since the matrix is diagonal the diagonal is equal to the eigenvalues of
  * the matrix.  Therefore you can show that the condition number measured in
  * any norm is max(diag(i))/min(diag(i)).  therefore we just need
  * to limit the smallest diagonal as diag(i) > max(diag(i)) / max_cond.
  */
class InitFinDiffReducedHessian_Step : public rSQPAlgo_Step {
public:

	///
	enum EInitializationMethod { SCALE_IDENTITY, SCALE_DIAGONAL
		, SCALE_DIAGONAL_ABS };

	///
	InitFinDiffReducedHessian_Step(
		  EInitializationMethod		initialization_method	= SCALE_IDENTITY
		, value_type				max_cond				= 1e+1
		, value_type				min_diag				= 1e-8
		, value_type				step_scale				= 1e-1			 );

	/// The initialization method for setting the diagonal
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EInitializationMethod,initialization_method)

	/// Maximum condition (l2 norm) for the intial matrix = (max(diag(i))/min(diag(i)).
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,max_cond)

	/// The absolute minimum value of a diagonal element
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,min_diag)

	/// The scaling of the step length u = step_scale / ||Z*e||inf
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,step_scale)

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	enum { NO_BASIS_UPDATED_YET = -1 };
	int num_basis_;
	quasi_newton_stats_iq_member	quasi_newton_stats_;

};	// end class ReducedHessianBFGS_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H