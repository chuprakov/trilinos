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

#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOpDecl.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpBaseDecl.hpp"
#include "Thyra_VectorDecl.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_LinearOperatorDecl.hpp"
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultBlockedLinearOpDecl.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Meros_LSCPreconditionerFactory.h"
#include "Meros_AztecSolveStrategy.hpp"
#include "Meros_InverseOperator.hpp"
#include "Meros_ZeroOperator.hpp"
#include "Meros_IdentityOperator.hpp"
#include "Meros_LinearSolver.hpp"
#include "Meros_PCDPreconditionerFactory.h"

#include "AztecOO.h"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"

using namespace Teuchos;
using namespace Thyra;
using namespace Meros;


// Constructors/initializers/accessors

// LSCPreconditionerFactory
// ::LSCPreconditionerFactory()
// {;}


LSCPreconditionerFactory
::LSCPreconditionerFactory(RefCountPtr<ParameterList> azFParams,
			   RefCountPtr<ParameterList> azBBtParams)
{
  azFParams_ = azFParams;
  azBBtParams_ = azBBtParams;
}


bool LSCPreconditionerFactory
::isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const
{
TEST_FOR_EXCEPT("LSCPreconditionerFactory::isCompatible is not implemented");
}


RefCountPtr<PreconditionerBase<double> > LSCPreconditionerFactory
::createPrec() const
{
  return rcp(new DefaultPreconditioner<double>());
}


void LSCPreconditionerFactory
::initializePrec(const RefCountPtr<const LinearOpSourceBase<double> > &opSrc,
		 PreconditionerBase<double> *prec,
		 const ESupportSolveUse supportSolveUse) const
{
  // Cast the general LinearOpSourceBase object to a LSCOperatorSource
  RefCountPtr<const LSCOperatorSource> lscOpSrcPtr 
    = rcp_dynamic_cast<const LSCOperatorSource>(opSrc);  
  
  // Retrieve operators from the LSC operator source
  ConstLinearOperator<double> blockOp = lscOpSrcPtr->getSaddleOp();
  ConstLinearOperator<double> F = blockOp.getBlock(0,0);
  ConstLinearOperator<double> Bt = blockOp.getBlock(0,1);
  ConstLinearOperator<double> B = blockOp.getBlock(1,0);
  // This version of LSC assumes a stable discretization. 
  // Ignoring C block.

  // Builde F inverse operator Finv
  LinearSolveStrategy<double> azF 
    = new AztecSolveStrategy(*(azFParams_.get()));
  ConstLinearOperator<double> Finv 
    = new InverseOperator<double>(F, azF);

  // Builde BBt inverse operator BBtinv
  ConstLinearOperator<double> BBt = B * Bt;
  LinearSolveStrategy<double> azBBt 
    = new AztecSolveStrategy(*(azBBtParams_.get()));
  //  ConstLinearOperator<double> BBtinv 
  //    = new InverseOperator<double>(makeEpetraOperator(BBt), azBBt);
  ConstLinearOperator<double> BBtinv 
    = new InverseOperator<double>(BBt, azBBt);

  // Build identity matrices on the velocity and pressure spaces 
  ConstLinearOperator<double> Ivel = new IdentityOperator<double>(Bt.range());
  ConstLinearOperator<double> Ipress = new IdentityOperator<double>(Bt.domain());

  // Build zero operators. Need one that is pressure x velocity and
  // one that is velocity x pressure
  ConstLinearOperator<double> Z = new ZeroOperator<double>(Bt.domain(), 
						      Bt.range());
  ConstLinearOperator<double> Zt = new ZeroOperator<double>(Bt.range(), 
						       Bt.domain());

  // Build Quinv
  // Setting Quinv = I for now. Will not be identity for some discretizations.
  ConstLinearOperator<double> Quinv = Ivel;

  ConstLinearOperator<double> BFBt = B * F * Bt;
  ConstLinearOperator<double> Xinv = BBtinv * BFBt *  BBtinv;

  // Build the 3 block operators for the preconditioner
  ConstLinearOperator<double> P1 = block2x2(Finv, Zt, Z, Ipress);
  ConstLinearOperator<double> P2 = block2x2(Ivel, (-1.0)*Bt, Z, Ipress);
  ConstLinearOperator<double> P3 = block2x2(Ivel, Zt, Z, (-1.0)*Xinv);

  //  ConstLinearOperator<double> LSCprec = makeEpetraOperator(P1 * P2 * P3);
  ConstLinearOperator<double> LSCprec = P1 * P2 * P3;

  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);

  (*defaultPrec).initializeRight(LSCprec.constPtr());

//   if(hasQp_)
//     {
//       RefCountPtr<const LinearOpBase<double> > tmpQpOp = lscOpSrcPtr->getQp();
//       ConstLinearOperator<double> QpOp = tmpQpOp;
//     }
//   else
//     {;
//     }

}


void LSCPreconditionerFactory
::uninitializePrec(PreconditionerBase<double> *prec,
		   RefCountPtr<const LinearOpSourceBase<double> > *fwdOp,
		   ESupportSolveUse *supportSolveUse) const
{
TEST_FOR_EXCEPT("LSCPreconditionerFactory::uninitializePrec not implemented");
}


// Overridden from ParameterListAcceptor

void LSCPreconditionerFactory
::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  // Don't know how to validate an ML list
  //  paramList->validateParameters(*this->getValidParameters(),1);
  paramList_ = paramList;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
LSCPreconditionerFactory::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
LSCPreconditionerFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
LSCPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
LSCPreconditionerFactory::getValidParameters() const
{
  // if(!validPL_.get()) {
//     validPL_ = defaultParameters(ML_DomainDecomposition);
//     // Todo: the above is not really all of the valid parameters.  We need to
//     // get ML to generate the entire list!
//   }
  return validPL_;
}
