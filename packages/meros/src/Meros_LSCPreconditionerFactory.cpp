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

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_InverseLinearOperator.hpp"

#include "Meros_LSCPreconditionerFactory.h"
#include "Meros_PCDPreconditionerFactory.h"

using namespace Meros;

// Constructors/initializers/accessors

// LSCPreconditionerFactory
// ::LSCPreconditionerFactory()
// {;}


LSCPreconditionerFactory
::LSCPreconditionerFactory(
  RCP<const LinearOpWithSolveFactoryBase<double> >   const&  FSolveStrategy,
  RCP<const LinearOpWithSolveFactoryBase<double> >   const&  BBtSolveStrategy
  )
  :FSolveStrategy_(FSolveStrategy),
   BBtSolveStrategy_(BBtSolveStrategy)
{}

bool LSCPreconditionerFactory
::isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const
{
  TEUCHOS_TEST_FOR_EXCEPT("LSCPreconditionerFactory::isCompatible is not implemented");
  return(true);
}

RCP<PreconditionerBase<double> > LSCPreconditionerFactory
::createPrec() const
{
  return rcp(new DefaultPreconditioner<double>());
}

void LSCPreconditionerFactory
::initializePrec(const RCP<const LinearOpSourceBase<double> > &opSrc,
		 PreconditionerBase<double> *prec,
		 const ESupportSolveUse supportSolveUse) const
{
  // Cast the general LinearOpSourceBase object to a LSCOperatorSource
  RCP<const LSCOperatorSource> lscOpSrcPtr 
    = rcp_dynamic_cast<const LSCOperatorSource>(opSrc);  
  
  // Retrieve operators from the LSC operator source
  ConstLinearOperator<double> blockOp = lscOpSrcPtr->getSaddleOp();
  ConstLinearOperator<double> F = blockOp.getBlock(0,0);
  ConstLinearOperator<double> Bt = blockOp.getBlock(0,1);
  ConstLinearOperator<double> B = blockOp.getBlock(1,0);
  // This version of LSC assumes a stable discretization. 
  // Ignoring C block.

  // Builde F inverse operator Finv
  ConstLinearOperator<double>
    Finv = inverse(*FSolveStrategy_,F,Thyra::IGNORE_SOLVE_FAILURE);

  // Builde BBt inverse operator BBtinv
  ConstLinearOperator<double> BBt = B * Bt;
  ConstLinearOperator<double>
    BBtinv = inverse(*BBtSolveStrategy_,BBt,Thyra::IGNORE_SOLVE_FAILURE);

  // Build identity matrices on the velocity and pressure spaces 
  ConstLinearOperator<double> Ivel = Thyra::identity<double>(Bt.range());
  ConstLinearOperator<double> Ipress = Thyra::identity<double>(Bt.domain());

  // Build zero operators. Need one that is pressure x velocity and
  // one that is velocity x pressure
  ConstLinearOperator<double> zero;

  // Build Quinv
  // Setting Quinv = I for now. Will not be identity for some discretizations.
  ConstLinearOperator<double> Quinv = Ivel;

  ConstLinearOperator<double> BFBt = B * F * Bt;
  ConstLinearOperator<double> Xinv = BBtinv * BFBt *  BBtinv;

  // Build the 3 block operators for the preconditioner
  ConstLinearOperator<double> P1 = block2x2(  Finv,    zero,
                                              zero,    Ipress        );
  ConstLinearOperator<double> P2 = block2x2(  Ivel,    (-1.0)*Bt,
                                              zero,    Ipress        );
  ConstLinearOperator<double> P3 = block2x2(  Ivel,    zero, 
                                              zero,    (-1.0)*Xinv   );

  //  ConstLinearOperator<double> LSCprec = makeEpetraOperator(P1 * P2 * P3);
  ConstLinearOperator<double> LSCprec = P1 * P2 * P3;

  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);

  (*defaultPrec).initializeRight(LSCprec.constPtr());

//   if(hasQp_)
//     {
//       RCP<const LinearOpBase<double> > tmpQpOp = lscOpSrcPtr->getQp();
//       ConstLinearOperator<double> QpOp = tmpQpOp;
//     }
//   else
//     {;
//     }

}


void LSCPreconditionerFactory
::uninitializePrec(PreconditionerBase<double> *prec,
		   RCP<const LinearOpSourceBase<double> > *fwdOp,
		   ESupportSolveUse *supportSolveUse) const
{
TEUCHOS_TEST_FOR_EXCEPT("LSCPreconditionerFactory::uninitializePrec not implemented");
}


// Overridden from ParameterListAcceptor

void LSCPreconditionerFactory
::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  // Don't know how to validate an ML list
  //  paramList->validateParameters(*this->getValidParameters(),1);
  paramList_ = paramList;
}

Teuchos::RCP<Teuchos::ParameterList>
LSCPreconditionerFactory::getNonconstParameterList()
{
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList>
LSCPreconditionerFactory::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList>
LSCPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList>
LSCPreconditionerFactory::getValidParameters() const
{
  // if(!validPL_.get()) {
//     validPL_ = defaultParameters(ML_DomainDecomposition);
//     // Todo: the above is not really all of the valid parameters.  We need to
//     // get ML to generate the entire list!
//   }
  return validPL_;
}
