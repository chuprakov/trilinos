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

#include "NSBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"

#include "RightBlockNSOperatorSource.h"
#include "KayLoghinRightOperatorSource.h"

using namespace SPP;
using namespace TSF;
using std::string;

NSBlockPreconditionerFactory
::NSBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac)
 	: TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac)
{}

TSFPreconditioner NSBlockPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& op) const
{
  TSFError::raise("NSBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFLinearOperator instead of TSFOperatorSource");
  return 0;
}

TSFPreconditioner NSBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  // Kay & Loghin style block preconditioner
  // 
  //      |  inv(F)   0   | |   I    -Bt  | |   I      0     |
  //      |    0      I   | |   0     I   | |   0   -inv(X)  |
  //

  // Get the saddle point operator S
  TSFLinearOperator S = opSource.getOp();
  
  // Get the blocks from S
  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Finv = F.inverse(Fsolver_);

  TSFLinearOperator Bt = S.getBlock(0,1);
  //	TSFLinearOperator C = S.getBlock(1,0);
  	
  // Build the Schur complement approx using the SchurFactory
  TSFLinearOperator Xinv = sfac_.getSchurInvApprox(opSource);

  // Make identity matrices on the right spaces
  TSFLinearOperator I00 = new TSFIdentityOperator(F.domain());
	TSFLinearOperator I11= new TSFIdentityOperator(Bt.domain());

  // Make block operators for the three preconditioner factors
	TSFLinearOperator P1 = new TSFBlockLinearOperator(S.domain(), S.range());
	TSFLinearOperator P2 = new TSFBlockLinearOperator(S.domain(), S.range());
	TSFLinearOperator P3 = new TSFBlockLinearOperator(S.domain(), S.range());

  // Put the pieces together into the composed block preconditioner

  P1.setBlock(0, 0, Finv);
	P1.setBlock(1, 1, I11);

	P2.setBlock(0, 0, I00);
	P2.setBlock(0, 1, -Bt);
	P2.setBlock(1, 1, I11);

	P3.setBlock(0, 0, I00);
	P3.setBlock(1, 1, -Xinv);

  return new GenericRightPreconditioner(P1*P2*P3);

}


string NSBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





