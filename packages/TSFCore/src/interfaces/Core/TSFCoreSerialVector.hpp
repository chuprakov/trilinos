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

// /////////////////////////////////////////////////////////////////
// TSFCoreSerialVector.hpp

#ifndef TSFCORE_VECTOR_SERIAL_HPP
#define TSFCORE_VECTOR_SERIAL_HPP

#include "TSFCoreSerialVectorDecl.hpp"
#include "TSFCoreSerialVectorBase.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVector<Scalar>::SerialVector(
	const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
	)
{
	this->initialize(vecSpc);
}

template<class Scalar>
SerialVector<Scalar>::SerialVector(
	const Index dim
	)
{
	this->initialize(dim);
}

template<class Scalar>
SerialVector<Scalar>::SerialVector(
	const Teuchos::RefCountPtr<Scalar>                      &v
	,const Index                                            vs
	,const Index                                            dim
	,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
	)
{
	this->initialize(v,vs,dim,vecSpc);
}

template<class Scalar>
void SerialVector<Scalar>::initialize(
	const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
	)
{
	const Index dim = vecSpc->dim();
	this->initialize(
		Teuchos::rcp( new Scalar[dim], Teuchos::DeallocArrayDelete<Scalar>(), true )
		,1
		,dim
		,vecSpc
		);
}

template<class Scalar>
void SerialVector<Scalar>::initialize(
	const Index dim
	)
{
	this->initialize(
		Teuchos::rcp( new Scalar[dim], Teuchos::DeallocArrayDelete<Scalar>(), true )
		,1
		,dim
		);
}

template<class Scalar>
void SerialVector<Scalar>::initialize(
	const Teuchos::RefCountPtr<Scalar>                      &v
	,const Index                                            vs
	,const Index                                            dim
	,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
	)
{
	if(vecSpc.get()) {
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vecSpc.get()!=NULL && dim != vecSpc->dim(), std::invalid_argument, "SerialVector<Scalar>::initialize(...): Error!" );
#endif
		space_serial_ = vecSpc;
	}
	else {
		space_serial_ = Teuchos::rcp(new SerialVectorSpace<Scalar>(dim));
	}
	v_       = v;
	vs_      = vs;
	dim_     = dim;
}

// Overridden from SerialVectorBase

template<class Scalar>
void SerialVector<Scalar>::getData( Scalar** values, Index* stride )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT(values==NULL || stride==NULL);
#endif
	*values = v_.get();
	*stride = vs_;
}

template<class Scalar>
void SerialVector<Scalar>::getData( const Scalar** values, Index* stride ) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT(values==NULL || stride==NULL);
#endif
	*values = v_.get();
	*stride = vs_;
}

// Overridden from Vector

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
SerialVector<Scalar>::space() const
{
	return space_serial_;
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_HPP
