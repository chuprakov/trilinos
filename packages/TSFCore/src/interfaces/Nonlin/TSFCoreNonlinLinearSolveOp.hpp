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
// TSFCoreNonlinLinearSolveOp.hpp

#ifndef TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP
#define TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP

#include "TSFCoreNonlinLinearSolveOpDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
ENonsingStatus LinearSolveOp<Scalar>::nonsingStatus() const
{
	return OP_SINGULARITY_UNKNOWN;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearSolveOp<Scalar> >
LinearSolveOp<Scalar>::clone_lso() const
{
	return Teuchos::null;
}

template<class Scalar>
void LinearSolveOp<Scalar>::solve(
	const ETransp                          M_trans
	,const MultiVector<Scalar>             &Y
	,MultiVector<Scalar>                   *X
	,const Scalar                          alpha
	,Solvers::ConvergenceTester<Scalar>    *convTester
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT(X==NULL);
	TSFCORE_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("Vector<Scalar>::apply()",*this,M_trans,*X,&Y);
#endif
	const VectorSpace<Scalar> &space_mv_rows = *Y.domain();
	const Index               num_mv_cols    = space_mv_rows.dim();
	//
	// Here we will solve the linear systems one at a time as:
	//
	//    op(M)*X(j) = alpha*Y(j))
	//
	bool scale_y = (alpha != 1.0);
	Teuchos::RefCountPtr<Vector<Scalar> > tmp = ( scale_y ? Y.range()->createMember() : Teuchos::null );
	for( Index j = 1; j <= num_mv_cols; ++j ) {
		// get tmp = alpha*Y.col(j) (but only scale if alpha != 1.0)
		if( scale_y ) {
			*tmp = *Y.col(j);
			Vt_S(tmp.get(),alpha);
		}
		else {
			tmp = Teuchos::rcp_const_cast<Vector<Scalar> >(Y.col(j));
		}
		// Solve op(M)*X(j) = alpha*Y(j)
		this->solve(M_trans,*tmp,X->col(j).get(),convTester);
	}
}

} // namespace Nonlin
} // namespace TSFCore

#endif /// TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP
