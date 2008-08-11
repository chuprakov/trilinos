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
#include "Meros_PCDOperatorSource.h"
#include "Meros_IdentityOperator.hpp"

using namespace Thyra;
using namespace Meros;

  
// Constructors/initializers/accessors
  

PCDOperatorSource::PCDOperatorSource()
{;}
  

PCDOperatorSource
::PCDOperatorSource(ConstLinearOperator<double> op,
		    ConstLinearOperator<double> Fp,
		    ConstLinearOperator<double> Ap)
{
  ConstLinearOperator<double> Bt = op.getBlock(0,1);
  ConstLinearOperator<double> Qp = new IdentityOperator<double>(Bt.domain());

  op_.initialize(op.constPtr());
  Fp_.initialize(Fp.constPtr());
  Ap_.initialize(Ap.constPtr());
  Qp_.initialize(Qp.constPtr());
}
  

PCDOperatorSource
::PCDOperatorSource(ConstLinearOperator<double> op,
		    ConstLinearOperator<double> Fp,
		    ConstLinearOperator<double> Ap,
		    ConstLinearOperator<double> Qp)
{
  op_.initialize(op.constPtr());
  Fp_.initialize(Fp.constPtr());
  Ap_.initialize(Ap.constPtr());
  Qp_.initialize(Qp.constPtr());
}
  

PCDOperatorSource
::PCDOperatorSource(Epetra_RowMatrix* S00epetra,
 		    Epetra_RowMatrix* S01epetra,
		    Epetra_RowMatrix* S10epetra,
		    Epetra_RowMatrix* S11epetra,
		    Epetra_RowMatrix* Fpepetra,
		    Epetra_RowMatrix* Apepetra)
{

  // Convert to LinearOperators
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

  // Build block matrix
  ConstLinearOperator<double> S = block2x2(S00, S01, S10,S11);

  // Convert Fp and Ap to Thyra LinearOperators
  RCP<const LinearOpBase<double> >
    tmpFp = Thyra::epetraLinearOp(rcp(Fpepetra,false));
  ConstLinearOperator<double> Fp = tmpFp;

  RCP<const LinearOpBase<double> >
    tmpAp = Thyra::epetraLinearOp(rcp(Apepetra,false));
  ConstLinearOperator<double> Ap = tmpAp;

  // Qp = I if not given in constructor
  ConstLinearOperator<double> Qp = new IdentityOperator<double>(S01.domain());

  op_.initialize(S.constPtr());
  Fp_.initialize(Fp.constPtr());
  Ap_.initialize(Ap.constPtr());
  Qp_.initialize(Qp.constPtr());
}
  

PCDOperatorSource
::PCDOperatorSource(Epetra_RowMatrix* S00epetra,
 		    Epetra_RowMatrix* S01epetra,
		    Epetra_RowMatrix* S10epetra,
		    Epetra_RowMatrix* S11epetra,
		    Epetra_RowMatrix* Fpepetra,
		    Epetra_RowMatrix* Apepetra,
		    Epetra_RowMatrix* Qpepetra)
{

  // Convert to Thyra LinearOperators
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

  // Build block operator
  ConstLinearOperator<double> S = block2x2(S00, S01, S10,S11);

  // Convert Fp, Ap, and Qp to Thyra LinearOperators
  RCP<const LinearOpBase<double> >
    tmpFp = Thyra::epetraLinearOp(rcp(Fpepetra,false));
  ConstLinearOperator<double> Fp = tmpFp;

  RCP<const LinearOpBase<double> >
    tmpAp = Thyra::epetraLinearOp(rcp(Apepetra,false));
  ConstLinearOperator<double> Ap = tmpAp;

  RCP<const LinearOpBase<double> >
    tmpQp = Thyra::epetraLinearOp(rcp(Qpepetra,false));
  ConstLinearOperator<double> Qp = tmpQp;

  op_.initialize(S.constPtr());
  Fp_.initialize(Fp.constPtr());
  Ap_.initialize(Ap.constPtr());
  Qp_.initialize(Qp.constPtr());
}
  



// void PCDOperatorSource
// ::initialize(ConstLinearOperator<double> op,
// 	     ConstLinearOperator<double> Fp,
// 	     ConstLinearOperator<double> Ap)
// {
//   op_.initialize(op.constPtr());
//   Fp_.initialize(Fp.constPtr());
//   Ap_.initialize(Ap.constPtr());
//   hasQp_ = false;
// }


void PCDOperatorSource
::initialize(ConstLinearOperator<double> op,
	     ConstLinearOperator<double> Fp,
	     ConstLinearOperator<double> Ap,
	     ConstLinearOperator<double> Qp)
{
  op_.initialize(op.constPtr());
  Fp_.initialize(Fp.constPtr());
  Ap_.initialize(Ap.constPtr());
  Qp_.initialize(Qp.constPtr());
}


void PCDOperatorSource::uninitialize()
{
  op_.uninitialize();
}

// Overridden from LinearOpSourceBase

bool PCDOperatorSource::isOpConst() const
{
  return op_.isConst();
}

RCP<LinearOpBase<double> > PCDOperatorSource::getNonconstOp() 
{
  return op_.getNonconstObj();
}

RCP<const LinearOpBase<double> > PCDOperatorSource::getOp() const
{
  return op_.getConstObj();
}

ConstLinearOperator<double> PCDOperatorSource::getSaddleOp() const
{
  return op_.getConstObj();
}

ConstLinearOperator<double> PCDOperatorSource::getFp() const
{
  return Fp_.getConstObj();
}

ConstLinearOperator<double> PCDOperatorSource::getAp() const
{
  return Ap_.getConstObj();
}

ConstLinearOperator<double> PCDOperatorSource::getQp() const
{
  return Qp_.getConstObj();
}

// See support/operator_solve/client_support/Thyra_DefaultLinearOpSource.hpp

