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

// ///////////////////////////////////////////////
// Epetra_DiagonalOperator.cpp

#include "Epetra_DiagonalOperator.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace Epetra {

// Constructors / initializers / accessors

DiagonalOperator::DiagonalOperator()
{}

DiagonalOperator::DiagonalOperator(
	const Teuchos::RefCountPtr<const Epetra_Map>      &map
	,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
	)
{
	initialize(map,diag);
}

void DiagonalOperator::initialize(
	const Teuchos::RefCountPtr<const Epetra_Map>      &map
	,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( map.get()==NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( diag.get()==NULL, std::invalid_argument, "Error!" );
	// ToDo: Validate that the maps are compatible!
#endif // _DEBUG
	map_ = map;
	diag_ = diag;
	UseTranspose_ = false;
}

void DiagonalOperator::uninitialize(
	Teuchos::RefCountPtr<const Epetra_Map>      *map
	,Teuchos::RefCountPtr<const Epetra_Vector>  *diag
	)
{
	if(map) *map = map_;
	if(diag) *diag = diag_;
	map_ = Teuchos::null;
	diag_ = Teuchos::null;
}

// Overridden from Epetra_Operator

int DiagonalOperator::SetUseTranspose(bool UseTranspose)
{
	UseTranspose_ = UseTranspose; // Totally ignored!
	return 0;
}

int DiagonalOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	Y.Multiply( 1.0, *diag_, X, 0.0 );  // Y = Diag * X (element-wise multiplication)
	return 0;
}

int DiagonalOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	Y.ReciprocalMultiply( 1.0, X, *diag_, 0.0 );  // Y = X / Diag (element-wise division)
	return 0;
}

double DiagonalOperator::NormInf() const
{
	double nrm[1];
	diag_->NormInf(nrm);
	return nrm[0];
}

const char* DiagonalOperator::Label() const
{
	return NULL;
}

bool DiagonalOperator::UseTranspose() const
{
	return UseTranspose_;
}

bool DiagonalOperator::HasNormInf() const
{
	return true;
}

const Epetra_Comm&
DiagonalOperator::Comm() const
{
	return map_->Comm();
}

const Epetra_Map&
DiagonalOperator::OperatorDomainMap() const
{
	return *map_;
}

const Epetra_Map&
DiagonalOperator::OperatorRangeMap() const
{
	return *map_;
}

} // namespace Epetra
