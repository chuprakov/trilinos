// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversIterativeLinearSolverDecl.hpp

#ifndef TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_DECL_HPP
#define TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_DECL_HPP

#include "TSFCoreSolversSolverState.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Strategy interface for preconditioned iterative linear solvers.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class IterativeLinearSolver : virtual public SolverState<Scalar> {
public:

	/** @name Pure virtual methods that must be overridden */
	//@{

 	///
	virtual bool adjointRequired() const = 0;

	///
	virtual SolveReturn solve(
		const LinearOp<Scalar>               &M
		,const ETransp                       M_trans
		,const MultiVector<Scalar>           &Y
		,MultiVector<Scalar>                 *X
		,const Scalar                        alpha                   = 1.0
		,const int                           max_iter                = DEFAULT_MAX_ITER
		,ConvergenceTester<Scalar>           *convTester             = NULL
		,const LinearOp<Scalar>              *M_tilde_left_inv       = NULL
		,const ETransp                       M_tilde_left_inv_trans  = NOTRANS
		,const LinearOp<Scalar>              *M_tilde_right_inv      = NULL
		,const ETransp                       M_tilde_right_inv_trans = NOTRANS
		) const = 0;

	//@}

	/** @name Virtual methods with default implementations */
	//@{

	///
	virtual SolveReturn solve(
		const LinearOp<Scalar>               &M
		,const ETransp                       M_trans
		,const Vector<Scalar>                &y
		,Vector<Scalar>                      *x
		,const int                           max_iter                = DEFAULT_MAX_ITER
		,ConvergenceTester<Scalar>           *convTester             = NULL
		,const LinearOp<Scalar>              *M_tilde_left_inv       = NULL
		,const ETransp                       M_tilde_left_inv_trans  = NOTRANS
		,const LinearOp<Scalar>              *M_tilde_right_inv      = NULL
		,const ETransp                       M_tilde_right_inv_trans = NOTRANS
		) const;

	///
	virtual Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> > clone() const;

	//@}

};	// end class IterativeLinearSolver

} // namespace Solvers
} // namespace TSFCore

#endif	// TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_DECL_HPP
