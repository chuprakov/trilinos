// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Meros_SIMPLEOperatorSource.h"
#include "Meros_IdentityOperator.hpp"

using namespace Thyra;
using namespace Meros;

  
// Constructors/initializers/accessors
  

SIMPLEOperatorSource::SIMPLEOperatorSource()
{;}
  

SIMPLEOperatorSource
::SIMPLEOperatorSource(ConstLinearOperator<double> op)
{
  ConstLinearOperator<double> F = op.getBlock(0,0);
  // ConstLinearOperator<double> Qu = new IdentityOperator<double>(F.domain());

  op_.initialize(op.constPtr());
  // Qu_.initialize(Qu.constPtr());
}
  

// SIMPLEOperatorSource
// ::SIMPLEOperatorSource(ConstLinearOperator<double> op,
// 		       ConstLinearOperator<double> Qu)
// {
//   op_.initialize(op.constPtr());
//   Qu_.initialize(Qu.constPtr());
// }
  

SIMPLEOperatorSource
::SIMPLEOperatorSource(Epetra_RowMatrix* S00epetra,
		       Epetra_RowMatrix* S01epetra,
		       Epetra_RowMatrix* S10epetra,
		       Epetra_RowMatrix* S11epetra)
{
  // convert to LinearOperators, build block matrix, and initialize
  RCP<const LinearOpBase<double> >
    tmpS00 = Thyra::epetraLinearOp(rcp(S00epetra,false));
  ConstLinearOperator<double> S00 = tmpS00;

  RCP<const LinearOpBase<double> >
    tmpS01 = Thyra::epetraLinearOp(rcp(S01epetra,false));
  ConstLinearOperator<double> S01 = tmpS01;

  RCP<const LinearOpBase<double> >
    tmpS10 = Thyra::epetraLinearOp(rcp(S10epetra,false));
  ConstLinearOperator<double> S10 = tmpS10;

  RCP<const LinearOpBase<double> >
    tmpS11 = Thyra::epetraLinearOp(rcp(S11epetra,false));
  ConstLinearOperator<double> S11 = tmpS11;

  ConstLinearOperator<double> S = block2x2(S00, S01, S10,S11);

  // ConstLinearOperator<double> Qu = new IdentityOperator<double>(S00.domain());

  op_.initialize(S.constPtr());
  // Qu_.initialize(Qu.constPtr());

}


// SIMPLEOperatorSource
// ::SIMPLEOperatorSource(Epetra_RowMatrix* S00epetra,
// 		       Epetra_RowMatrix* S01epetra,
// 		       Epetra_RowMatrix* S10epetra,
// 		       Epetra_RowMatrix* S11epetra,
// 		       Epetra_RowMatrix* Quepetra)
// {
//   // convert to LinearOperators, build block matrix, and initialize
//   RCP<const LinearOpBase<double> >
//     tmpS00 = rcp(new EpetraLinearOp(rcp(S00epetra,false)));
//   ConstLinearOperator<double> S00 = tmpS00;

//   RCP<const LinearOpBase<double> >
//     tmpS01 = rcp(new EpetraLinearOp(rcp(S01epetra,false)));
//   ConstLinearOperator<double> S01 = tmpS01;

//   RCP<const LinearOpBase<double> >
//     tmpS10 = rcp(new EpetraLinearOp(rcp(S10epetra,false)));
//   ConstLinearOperator<double> S10 = tmpS10;

//   RCP<const LinearOpBase<double> >
//     tmpS11 = rcp(new EpetraLinearOp(rcp(S11epetra,false)));
//   ConstLinearOperator<double> S11 = tmpS11;

//   RCP<const LinearOpBase<double> >
//     tmpQu = rcp(new EpetraLinearOp(rcp(Quepetra,false)));
//   ConstLinearOperator<double> Qu = tmpQu;

//   ConstLinearOperator<double> S = block2x2(S00, S01, S10,S11);

//   op_.initialize(S.constPtr());
//   Qu_.initialize(Qu.constPtr());

// }


// void SIMPLEOperatorSource
// ::initialize(ConstLinearOperator<double> op)
// {
//   op_.initialize(op.constPtr());
// }


// void SIMPLEOperatorSource
// ::initialize(ConstLinearOperator<double> op,
// 	     ConstLinearOperator<double> Qu)
// {
//   op_.initialize(op.constPtr());
//   Qu_.initialize(Qu.constPtr());
// }


void SIMPLEOperatorSource::uninitialize()
{
  op_.uninitialize();
}

// Overridden from LinearOpSourceBase

bool SIMPLEOperatorSource::isOpConst() const
{
  return op_.isConst();
}

RCP<LinearOpBase<double> > SIMPLEOperatorSource::getNonconstOp() 
{
  return op_.getNonconstObj();
}

RCP<const LinearOpBase<double> > SIMPLEOperatorSource::getOp() const
{
  return op_.getConstObj();
}

ConstLinearOperator<double> SIMPLEOperatorSource::getDinvOp() const
{
  return dinv_.getConstObj();
}

//ConstLinearOperator<double> SIMPLEOperatorSource::getQu() const
//{
//  return Qu_.getConstObj();
//}

// See support/operator_solve/client_support/Thyra_DefaultLinearOpSource.hpp

