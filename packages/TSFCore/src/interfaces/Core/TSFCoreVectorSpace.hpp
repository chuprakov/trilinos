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

// //////////////////////////////////////////////////////////////////////
// TSFCoreVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_HPP
#define TSFCORE_VECTOR_SPACE_HPP

#include "TSFCoreVectorSpaceDecl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreSerialVectorSpaceFactory.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreMultiVectorCols.hpp"

namespace TSFCore {

// Virtual functions with default implementations

template<class Scalar>
bool VectorSpace<Scalar>::isInCore() const
{
	return false;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceFactory<Scalar> >
VectorSpace<Scalar>::smallVecSpcFcty() const
{
	return Teuchos::rcp(new SerialVectorSpaceFactory<Scalar>());
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> > 
VectorSpace<Scalar>::createMembers(int numMembers) const
{
	return Teuchos::rcp(new MultiVectorCols<Scalar> (Teuchos::rcp(this,false),this->smallVecSpcFcty()->createVecSpc(numMembers)));
}

template<class Scalar>
Scalar VectorSpace<Scalar>::scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const
{
	const MultiVectorCols<Scalar>
		X( Teuchos::rcp( const_cast<Vector<Scalar>*>(&x), false ) ),
		Y( Teuchos::rcp( const_cast<Vector<Scalar>*>(&y), false ) );
	Scalar scalar_prods[1];
	this->scalarProds(X,Y,scalar_prods);
	return scalar_prods[0];
}

template<class Scalar>
void VectorSpace<Scalar>::scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const
{
	dot(X,Y,scalar_prods);
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
VectorSpace<Scalar>::clone() const
{
	return Teuchos::null;
}


} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_HPP
