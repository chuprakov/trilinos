// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversBiCGSolverDecl.hpp

#ifndef TSFCORE_SOLVERS_BICG_SOLVER_DECL_HPP
#define TSFCORE_SOLVERS_BICG_SOLVER_DECL_HPP

#include "TSFCoreSolversIterativeLinearSolver.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Implementation of a block preconditioned BiCG iterative solver.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class BiCGSolver : virtual public IterativeLinearSolver<Scalar> {
public:

	///
	/** Set the output steam for iterative algorithm.
	 */
	STANDARD_COMPOSITION_MEMBERS(std::ostream,out);

	///
	/** Set the default maximum number BiCG iterations to take
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );

	///
	/** Set the default maximum number BiCG iterations to take
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, default_max_iter );
	
	///
	/** Set the default solution tolerance to use
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, default_tol );

	///
	BiCGSolver(
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

	mutable Index                                    totalNumSystems_;
	mutable Index                                    currNumSystems_;
	mutable int                                      currIteration_;
	mutable std::vector<Index>                       activeSystems_;
	mutable Teuchos::RefCountPtr<const Solvers::Norm<Scalar> >
	    norm_;
	mutable Teuchos::RefCountPtr<MultiVector<Scalar> >
	    X_curr_, R_, R_tilde_, Q_, Q_tilde_, Z_, Z_tilde_, P_, P_tilde_;
	mutable std::vector<Scalar>
	    rho_, rho_old_, beta_, gamma_, alpha_;
	mutable bool norms_updated_;
	mutable std::vector<Scalar>
	    rel_err_denom_, norms_;

	// ///////////////////////////////
	// Private member functions
	
	///
	void doIteration(
		const LinearOp<Scalar>              &M
		,ETransp                            opM_notrans
		,ETransp                            opM_trans
		,MultiVector<Scalar>                *X
		,Scalar                             a
		,const LinearOp<Scalar>             *M_tilde_inv
		,ETransp                            opM_tilde_inv_notrans
		,ETransp                            opM_tilde_inv_trans
		) const;
	
	///
	void compress( bool isConverged[], MultiVector<Scalar>* X ) const;

};	// end class BiCGSolver

} // namespace Solvers
} // namespace TSFCore

#endif	// TSFCORE_SOLVERS_BICG_SOLVER_DECL_HPP
