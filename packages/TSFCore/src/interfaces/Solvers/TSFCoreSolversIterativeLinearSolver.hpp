// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversIterativeLinearSolver.hpp

#ifndef TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP
#define TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP

#include "TSFCoreSolversIterativeLinearSolverDecl.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
SolveReturn IterativeLinearSolver<Scalar>::solve(
	const LinearOp<Scalar>               &M
	,const ETransp                       M_trans
	,const Vector<Scalar>                &y
	,Vector<Scalar>                      *x
	,const int                           max_iter
	,ConvergenceTester<Scalar>           *convTester
	,const LinearOp<Scalar>              *M_tilde_left_inv
	,const ETransp                       M_tilde_left_inv_trans
	,const LinearOp<Scalar>              *M_tilde_right_inv
	,const ETransp                       M_tilde_right_inv_trans
	) const
{
	namespace mmp = MemMngPack;
	const MultiVectorCols<Scalar>  Y(mmp::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(mmp::rcp(x,false));
	return solve(
		M,M_trans,Y,&X,1.0,max_iter,convTester
		,M_tilde_left_inv,M_tilde_left_inv_trans
		,M_tilde_right_inv,M_tilde_right_inv_trans
		);
}

template<class Scalar>
MemMngPack::ref_count_ptr<const IterativeLinearSolver<Scalar> >
IterativeLinearSolver<Scalar>::clone() const
{
	return MemMngPack::null;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP
