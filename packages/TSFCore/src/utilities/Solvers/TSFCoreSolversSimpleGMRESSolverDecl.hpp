// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversSimpleGMRESSolverDecl.hpp

#ifndef TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_DECL_HPP
#define TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_DECL_HPP

#include "TSFCoreSolversIterativeLinearSolver.hpp"
#include "TSFCoreSolversGmresSolver.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Implementation of a simple GMRES iterative solver.
 *
 * ToDo: Finish documentation!
 *
 * Note, This implementation currently ignores the convergence tester!
 */
template<class Scalar>
class SimpleGMRESSolver : virtual public IterativeLinearSolver<Scalar> {
public:

	///
	/** Set the output steam for iterative algorithm.
	 */
	STANDARD_COMPOSITION_MEMBERS(std::ostream,out);

	///
	/** Set the default maximum number GMRES iterations to take
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );

	///
	/** Set the default maximum number GMRES iterations to take
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, default_max_iter );
	
	///
	/** Set the default solution tolerance to use
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, default_tol );

	///
	SimpleGMRESSolver(
		const out_ptr_t   &out               = Teuchos::null
		,bool             dump_all           = false
		,int              default_max_iter   = 10000
		,Scalar           default_tol        = 1e-10
		);

	/** @name Overridden from SolverState (only to be called by ConvergenceTester objects) */
	//@{

	///
	Index totalNumSystems() const;
	///
	Index currNumSystems() const;
	///
	int currIteration() const;
	///
	void currActiveSystems( Index activeSystems[] ) const;
	///
	void currEstRelResidualNorms( Scalar norms[] ) const;

	//@}

	/** @name Overridden from IterativeLinearSolver */
	//@{

 	///
	bool adjointRequired() const;
	///
	SolveReturn solve(
		const LinearOp<Scalar>               &M
		,const ETransp                       M_trans
		,const MultiVector<Scalar>           &Y
		,MultiVector<Scalar>                 *X
		,const Scalar                        alpha
		,const int                           max_iter
		,ConvergenceTester<Scalar>           *convTester
		,const LinearOp<Scalar>              *M_tilde_left_inv
		,const ETransp                       M_tilde_left_inv_trans
		,const LinearOp<Scalar>              *M_tilde_right_inv
		,const ETransp                       M_tilde_right_inv_trans
		) const;
	///
	Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> > clone() const;

	//@}

private:

	// ///////////////////////////////
	// Private data members

	mutable GMRESSolver<Scalar>  solver_;

};	// end class SimpleGMRESSolver

} // namespace Solvers
} // namespace TSFCore

#endif	// TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_DECL_HPP
