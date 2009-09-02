// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_SINGLE_RHS_EUCLIDEAN_LINEAR_OP_BASE_HPP
#define THYRA_SINGLE_RHS_EUCLIDEAN_LINEAR_OP_BASE_HPP

#include "Thyra_SingleRhsEuclideanLinearOpBaseDecl.hpp"
#include "Thyra_SingleScalarEuclideanLinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_AssertOp.hpp"

namespace Thyra {

template<class Scalar>
void SingleRhsEuclideanLinearOpBase<Scalar>::euclideanApply(
  const EOpTransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  const VectorSpaceBase<Scalar> &space_mv_rows = *Y->domain();
  const Index num_mv_cols    = space_mv_rows.dim();
  for( Index j = 0; j < num_mv_cols; ++j )
    this->euclideanApply(M_trans,*X.col(j),Y->col(j).get(),alpha,beta);
}

}	// end namespace Thyra

#endif // THYRA_SINGLE_RHS_EUCLIDEAN_LINEAR_OP_BASE_HPP
