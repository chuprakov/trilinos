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
// TSFCoreEpetraVector.cpp

#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreMPIVectorBase.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreMPIVectorBase.hpp"
#include "Teuchos_TestForException.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace TSFCore {

// Constructors/initializers

EpetraVector::EpetraVector()
{}

EpetraVector::EpetraVector(
	const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
	)
{
	initialize(epetra_vec,epetra_vec_spc);
}

void EpetraVector::initialize(
	const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
	)
{
  using Teuchos::dyn_cast;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_vec.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
	TEST_FOR_EXCEPTION( !epetra_vec_spc.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
#endif
	epetra_vec_ = epetra_vec;
  epetra_vec_spc_ = epetra_vec_spc;
	updateMpiSpace();
}

void EpetraVector::setUninitialized(
	Teuchos::RefCountPtr<Epetra_Vector>              *epetra_vec
	,Teuchos::RefCountPtr<const EpetraVectorSpace>   *epetra_vec_spc
	)
{
	if(epetra_vec) *epetra_vec = epetra_vec_;
	if(epetra_vec_spc) *epetra_vec_spc = epetra_vec_spc_;

	epetra_vec_ = Teuchos::null;
	epetra_vec_spc_ = Teuchos::null;

	updateMpiSpace();
}

// Overridden from MPIVectorBase

Teuchos::RefCountPtr<const MPIVectorSpaceBase<EpetraVector::Scalar> >
EpetraVector::mpiSpace() const
{
	return epetra_vec_spc_; // will be null if this is uninitialized!
}

void EpetraVector::getLocalData( Scalar** localValues, Index* stride )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues==NULL || stride==NULL );
#endif
	*localValues = &(*epetra_vec_)[0];
	*stride = 1;
}

void EpetraVector::commitLocalData( Scalar* localValues )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues != &(*epetra_vec_)[0] );
#endif
	// Nothing to commit!
}

void EpetraVector::getLocalData( const Scalar** localValues, Index* stride ) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues==NULL || stride==NULL );
#endif
	*localValues = &(*epetra_vec_)[0];
	*stride = 1;
}

void EpetraVector::freeLocalData( const Scalar* localValues ) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues != &(*epetra_vec_)[0] );
#endif
	// Nothing to commit!
}

} // end namespace TSFCore
