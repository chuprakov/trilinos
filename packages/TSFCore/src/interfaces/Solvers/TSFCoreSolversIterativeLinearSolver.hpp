// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

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
#ifdef TSFCORE_VECTOR_DERIVE_FROM_MULTI_VECTOR
	return solve(
		M,M_trans,static_cast<const MultiVector<Scalar>&>(y),static_cast<MultiVector<Scalar>*>(x)
    ,1.0,max_iter,convTester
		,M_tilde_left_inv,M_tilde_left_inv_trans
		,M_tilde_right_inv,M_tilde_right_inv_trans
		);
#else
	const MultiVectorCols<Scalar>  Y(Teuchos::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(Teuchos::rcp(x,false));
	return solve(
		M,M_trans,Y,&X,1.0,max_iter,convTester
		,M_tilde_left_inv,M_tilde_left_inv_trans
		,M_tilde_right_inv,M_tilde_right_inv_trans
		);
#endif
}

template<class Scalar>
Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> >
IterativeLinearSolver<Scalar>::clone() const
{
	return Teuchos::null;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP
