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
// TSFCoreEpetraVectorSpaceFactory.cpp

#include "TSFCoreEpetraVectorSpaceFactory.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

namespace TSFCore {

// Constructors / initializers

EpetraVectorSpaceFactory::EpetraVectorSpaceFactory()
{}

EpetraVectorSpaceFactory::EpetraVectorSpaceFactory(
	const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
	)
{
	initialize(epetra_comm);
}

void EpetraVectorSpaceFactory::initialize(
	const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
	)
{
	epetra_comm_ = epetra_comm;
}

void EpetraVectorSpaceFactory::setUninitialized(
	Teuchos::RefCountPtr<const Epetra_Comm> *epetra_comm
	)
{
	if(epetra_comm) *epetra_comm = epetra_comm_;
	epetra_comm_ = Teuchos::null;
}

// Overridden from VectorSpaceFactory

Teuchos::RefCountPtr< const VectorSpace<EpetraVectorSpaceFactory::Scalar> >
EpetraVectorSpaceFactory::createVecSpc(int dim) const
{
	Teuchos::RefCountPtr<const Epetra_LocalMap>
		epetra_map = Teuchos::rcp(new Epetra_LocalMap(dim,0,*epetra_comm_));
	Teuchos::set_extra_data( epetra_comm_, "epetra_comm", &epetra_map );
	return Teuchos::rcp(new EpetraVectorSpace(epetra_map));
}

} // end namespace TSFCore
