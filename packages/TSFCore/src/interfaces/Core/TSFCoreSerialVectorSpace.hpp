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

// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_SERIAL_HPP
#define TSFCORE_VECTOR_SPACE_SERIAL_HPP

#include "TSFCoreSerialVectorSpaceDecl.hpp"
#include "TSFCoreSerialVectorSpaceBase.hpp"
#include "TSFCoreSerialVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVectorSpace<Scalar>::SerialVectorSpace( int dim )
{
	initialize(dim);
}

template<class Scalar>
void SerialVectorSpace<Scalar>::initialize( int dim )
{
	dim_ = dim;
}

// Overridden from VectorSpace

template<class Scalar>
Index SerialVectorSpace<Scalar>::dim() const
{
	return dim_;
}

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
SerialVectorSpace<Scalar>::createMember() const
{
	return Teuchos::rcp(new SerialVector<Scalar>(dim_));
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
SerialVectorSpace<Scalar>::clone() const
{
	return Teuchos::rcp(new SerialVectorSpace<Scalar>(*this));
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_SERIAL_HPP
