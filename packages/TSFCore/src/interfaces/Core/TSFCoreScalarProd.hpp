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

// ///////////////////////////////////////////////////////////////
// TSFCoreScalarProd.hpp

#ifndef TSFCORE_SCALAR_PROD_HPP
#define TSFCORE_SCALAR_PROD_HPP

#include "TSFCoreScalarProdDecl.hpp"
#include "TSFCoreAssertOp.hpp"
#include "TSFCoreLinearOp.hpp"

namespace TSFCore {

template<class Scalar>
Scalar ScalarProd<Scalar>::scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const
{
	Scalar scalar_prods[1];
#ifdef TSFCORE_VECTOR_DERIVE_FROM_MULTI_VECTOR
	this->scalarProds(
    static_cast<const MultiVector<Scalar>&>(x)
    ,static_cast<const MultiVector<Scalar>&>(y)
    ,scalar_prods
    );
#else
	const MultiVectorCols<Scalar>
		X( Teuchos::rcp( const_cast<Vector<Scalar>*>(&x), false ) ),
		Y( Teuchos::rcp( const_cast<Vector<Scalar>*>(&y), false ) );
	this->scalarProds(X,Y,scalar_prods);
#endif
	return scalar_prods[0];
}

template<class Scalar>
void ScalarProd<Scalar>::apply(
	const EuclideanLinearOpBase<Scalar>   &M
	,const ETransp                        M_trans
	,const MultiVector<Scalar>            &X
	,MultiVector<Scalar>                  *Y
	,const Scalar                         alpha
	,const Scalar                         beta
	) const
{
#ifdef _DEBUG
	TSFCORE_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("EuclideanLinearOpBase<Scalar>::euclideanApply(...)",M,M_trans,X,Y);
#endif
	const Index numMv = X.domain()->dim();
	for( int j = 1; j <= numMv; ++j )
		this->apply( M, M_trans, *X.col(j), &*Y->col(j), alpha, beta );
}

} // end namespace TSFCore

#endif  // TSFCORE_SCALAR_PROD_HPP
