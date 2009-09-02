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

//PCDPreconditionerFactory
//::PCDPreconditionerFactory()
//{;}

PCDPreconditionerFactory
::PCDPreconditionerFactory(
  RCP<const LinearOpWithSolveFactoryBase<double> > 
  const&  FSolveStrategy,
  RCP<const LinearOpWithSolveFactoryBase<double> >  
  const&  ApSolveStrategy
  )
  :FSolveStrategy_(FSolveStrategy),
   ApSolveStrategy_(ApSolveStrategy),
   QpSolveStrategy_(Teuchos::null)
{}


PCDPreconditionerFactory
::PCDPreconditionerFactory(
  RCP<const LinearOpWithSolveFactoryBase<double> >   
  const&  FSolveStrategy,
  RCP<const LinearOpWithSolveFactoryBase<double> >  
  const&  ApSolveStrategy,
  RCP<const LinearOpWithSolveFactoryBase<double> >  
  const&  QpSolveStrategy
  )
  :FSolveStrategy_(FSolveStrategy),
   ApSolveStrategy_(ApSolveStrategy),
   QpSolveStrategy_(QpSolveStrategy)
{}


// PCDPreconditionerFactory
// ::PCDPreconditionerFactory(RCP<ParameterList> azFParams,
// 			   RCP<ParameterList> azApParams,
// 			   RCP<ParameterList> azQpParams)
// {
//   azFParams_ = azFParams;
//   azApParams_ = azApParams;
//   azQpParams_ = azQpParams;
// }



bool PCDPreconditionerFactory
::isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const
{
  TEST_FOR_EXCEPT("PCDPreconditionerFactory::isCompatible is not implemented");
  return(true);
}


RCP<PreconditionerBase<double> > PCDPreconditionerFactory
::createPrec() const
{
  return rcp(new DefaultPreconditioner<double>());
}


void PCDPreconditionerFactory
::initializePrec(const RCP<const LinearOpSourceBase<double> > &opSrc,
		 PreconditionerBase<double> *prec,
		 const ESupportSolveUse supportSolveUse) const
{
  // Cast the LinearOpSourceBase object back to a PCDOperatorSource
  RCP<const PCDOperatorSource> pcdOpSrcPtr 
    = rcp_dynamic_cast<const PCDOperatorSource>(opSrc);

  // Retrieve block operator and subblocks from the PCD operator source
  //RCP<const LinearOpBase<double> > tmpBlockOp = pcdOpSrcPtr->getOp();
  ConstLinearOperator<double> blockOp = pcdOpSrcPtr->getSaddleOp();
  ConstLinearOperator<double> F = blockOp.getBlock(0,0);
  ConstLinearOperator<double> Bt = blockOp.getBlock(0,1);
  ConstLinearOperator<double> B = blockOp.getBlock(1,0);
  // don't need C block for PCD

  // Retrieve Fp, Ap, and Qp operators from PCD operator source.
  // Note that Qp is set to the Identity if not included in the 
  // operator source.
  ConstLinearOperator<double> Fp = pcdOpSrcPtr->getFp();
  ConstLinearOperator<double> Ap = pcdOpSrcPtr->getAp();
  ConstLinearOperator<double> Qp = pcdOpSrcPtr->getQp();
  
  // Builde F inverse operator Finv using the AztecOO solver param list.
  // LinearSolveStrategy<double> azF 
  //   = new AztecSolveStrategy(*(azFParams_.get()));
  ConstLinearOperator<double> 
    Finv = inverse(*FSolveStrategy_,F,Thyra::IGNORE_SOLVE_FAILURE);

  // Build Ap inverse operator Apinv using the AztecOO solver param list.
  //  LinearSolveStrategy<double> azAp 
  //   = new AztecSolveStrategy(*(azApParams_.get()));
  ConstLinearOperator<double> 
    Apinv = inverse(*ApSolveStrategy_,Ap,Thyra::IGNORE_SOLVE_FAILURE);

  // Build Qp inverse operator Qpinv using the AztecOO solver param list.
  // If no Qpinv solver was given, we'll use the Ap solver parameters.
  // LinearSolveStrategy<double> azQp;
  ConstLinearOperator<double> Qpinv; 
  if(QpSolveStrategy_.get() != NULL)
    {
      // use given Qp params
      // azQp = new AztecSolveStrategy(*(azQpParams_.get()));
      Qpinv = inverse(*QpSolveStrategy_,Qp,Thyra::IGNORE_SOLVE_FAILURE);
    }
  else
    {
      // use the Ap params
      // azQp = new AztecSolveStrategy(*(azApParams_.get()));
      Qpinv = inverse(*ApSolveStrategy_,Qp,Thyra::IGNORE_SOLVE_FAILURE);
    }

  // Build identity matrices on the velocity and pressure spaces 
  ConstLinearOperator<double> Ivel = Thyra::identity<double>(Bt.range());
  ConstLinearOperator<double> Ipress = Thyra::identity<double>(Bt.domain());

  // Build zero operators. Need one that is pressure x velocity and
  // one that is velocity x pressure
  ConstLinearOperator<double> zero;

  // Build the composed Schur complement approximation inverse
  ConstLinearOperator<double> Xinv = Qpinv * Fp * Apinv;

  // Build the 3 block operators for the preconditioner
  ConstLinearOperator<double> P1 = block2x2( Finv, zero, zero, Ipress );
  ConstLinearOperator<double> P2 = block2x2( Ivel, (-1.0)*Bt, zero, Ipress  );
  ConstLinearOperator<double> P3 = block2x2( Ivel,   zero, zero, (-1.0)*Xinv );

  ConstLinearOperator<double> PCDprec = P1 * P2 * P3;

  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);

  (*defaultPrec).initializeRight(PCDprec.constPtr());

}


void PCDPreconditionerFactory
::uninitializePrec(PreconditionerBase<double> *prec,
		   RCP<const LinearOpSourceBase<double> > *fwdOp,
		   ESupportSolveUse *supportSolveUse) const
{
TEST_FOR_EXCEPT("PCDPreconditionerFactory::uninitializePrec not implemented");
}


// Overridden from ParameterListAcceptor

void PCDPreconditionerFactory
::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  // Don't really have meros parameter lists yet.
  //  paramList->validateParameters(*this->getValidParameters(),1);
  paramList_ = paramList;
}

Teuchos::RCP<Teuchos::ParameterList>
PCDPreconditionerFactory::getNonconstParameterList()
{
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList>
PCDPreconditionerFactory::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList>
PCDPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList>
PCDPreconditionerFactory::getValidParameters() const
{
  // if(!validPL_.get()) { validPL_ =
  // defaultParameters(ML_DomainDecomposition); // Todo: the above is
  // not really all of the valid parameters.  We need to // get ML to
  // generate the entire list!  }
  return validPL_;
}
