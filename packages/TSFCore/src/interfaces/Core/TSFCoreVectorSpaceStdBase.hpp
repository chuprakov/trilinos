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
// TSFCoreVectorSpaceStdBase.hpp

#ifndef TSFCORE_VECTOR_SPACE_STD_BASE_HPP
#define TSFCORE_VECTOR_SPACE_STD_BASE_HPP

#include "TSFCoreVectorSpaceStdBaseDecl.hpp"
#include "TSFCoreDotProd.hpp"
#include "TSFCoreAssertOp.hpp"

namespace TSFCore {

// Constructors / initializers

template<class Scalar>
VectorSpaceStdBase<Scalar>::VectorSpaceStdBase()
	:scalarProd_(Teuchos::rcp(new DotProd<Scalar>()))
{}
	
template<class Scalar>
VectorSpaceStdBase<Scalar>::VectorSpaceStdBase( const Teuchos::RefCountPtr<const ScalarProd<Scalar> > &scalarProd )
	:scalarProd_(scalarProd)
{
	TEST_FOR_EXCEPT( scalarProd.get()==NULL );
}

template<class Scalar>
void VectorSpaceStdBase<Scalar>::setScalarProd( const Teuchos::RefCountPtr<const ScalarProd<Scalar> > &scalarProd )
{
	TEST_FOR_EXCEPT( scalarProd.get()==NULL );
	scalarProd_ = scalarProd;
}

template<class Scalar>
Teuchos::RefCountPtr<const ScalarProd<Scalar> >
VectorSpaceStdBase<Scalar>::getScalarProd() const
{
	return scalarProd_;
}

// Overridden from VectorSpace

template<class Scalar>
Scalar VectorSpaceStdBase<Scalar>::scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const
{
#ifdef _DEBUG
	TSFCORE_ASSERT_VEC_SPACES("VectorSpaceStdBase<Scalar>::scalarProd(...)",*x.space(),*this);
	TSFCORE_ASSERT_VEC_SPACES("VectorSpaceStdBase<Scalar>::scalarProd(...)",*y.space(),*this);
#endif
	return scalarProd_->scalarProd(x,y);
}

template<class Scalar>
void VectorSpaceStdBase<Scalar>::scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const
{
#ifdef _DEBUG
	TSFCORE_ASSERT_VEC_SPACES("VectorSpaceStdBase<Scalar>::scalarProds(...)",*X.range(),*this);
	TSFCORE_ASSERT_VEC_SPACES("VectorSpaceStdBase<Scalar>::scalarProds(...)",*Y.range(),*this);
	TSFCORE_ASSERT_VEC_SPACES("VectorSpaceStdBase<Scalar>::scalarProds(...)",*X.domain(),*Y.domain());
#endif
	scalarProd_->scalarProds(X,Y,scalar_prods);
}

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_STD_BASE_HPP
