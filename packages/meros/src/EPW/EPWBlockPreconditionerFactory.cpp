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

#include "EPWBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"
#include "RightBlockNSOperatorSource.h"
#include "BlockForwardsolver.h"
#include "TSFZeroOperator.h"
#include "Aztec2TSF.h"
#include "Epetra_CrsMatrix.h"
#include "PetraMatrix.h"


using namespace Meros;
using namespace TSF;
using std::string;

//Need to update to reflect the Schur solver (pick a default solver)
EPWBlockPreconditionerFactory
::EPWBlockPreconditionerFactory(const TSFLinearSolver& blockSolver)
 	: TSFPreconditionerFactoryBase(), 
	  blockSolver_(blockSolver),
	  SchurSolver_(blockSolver)
{}

//Need to update to reflect the Schur solver (pick a default solver)
EPWBlockPreconditionerFactory
::EPWBlockPreconditionerFactory(const TSFLinearSolver& blockSolver,
				const TSFLinearSolver& SchurSolver)
 	: TSFPreconditionerFactoryBase(), 
	  blockSolver_(blockSolver),
	  SchurSolver_(SchurSolver)
{}


TSFPreconditioner EPWBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  TSFError::raise("EPWBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFOperatorSource instead of TSFLinearOperator");
  return 0;
}

TSFPreconditioner EPWBlockPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& op) const
{
  // block preconditioner
  // no Schur approximation yet
  // 
  //      |  inv(MpN)   0   | |   I    -N   | |   I      0     |
  //      |    0      I     | |   0     I   | |   0   -inv(X)  |
  // X = S = Schur complement

  // Get the saddle point operator S
  //  int numMyRows = S_epet->NumMyRows();
  //  cerr << "numMyRows = " << numMyRows << endl;
  cerr << "problem size = " << op.domain().dim() << endl;

  //   Epetra_Map map (numGlobalElements, 0, comm);
  //  Epetra_CrsMatrix A (Copy, map, 
  // GetCrsMatrix("matrixA.data", testA, 0);


  // MpN = M+N, but we don't have M separately
  TSFLinearOperator MpN = op.getBlock(0,0);
  TSFLinearOperator Nt = op.getBlock(0,1);
  // These should be the same as MpN and N (except for BCs)
  TSFLinearOperator N = op.getBlock(1,0);
  TSFLinearOperator MpN2 = op.getBlock(1,1);

  TSFLinearOperator MpNinv = MpN.inverse(blockSolver_);

  TSFLinearOperator X = MpN2 - N*MpNinv*Nt;
  TSFLinearOperator Xinv = X.inverse(SchurSolver_);

  // Make identity matrices on the right spaces
  TSFLinearOperator Ix = new TSFIdentityOperator(MpN.domain());
  TSFLinearOperator Iy = new TSFIdentityOperator(Nt.domain());
  
  // Make block operators for the three preconditioner factors
  TSFLinearOperator P1 = new TSFBlockLinearOperator(op.domain(), op.range());
  TSFLinearOperator P2 = new TSFBlockLinearOperator(op.domain(), op.range());
  TSFLinearOperator P3 = new TSFBlockLinearOperator(op.domain(), op.range());

  // Put the pieces together into the composed block preconditioner
  P1.setBlock(0,0,MpNinv);
  P1.setBlock(1,1,Iy);
  
  P2.setBlock(0,0,Ix);
  P2.setBlock(0,1,-Nt);
  P2.setBlock(1,1,Iy);
  
  P3.setBlock(0,0,Ix);
  P3.setBlock(1,1,Xinv);

// Testing just to see how we do on the (0,0) block
//  P1.setBlock(0,0,MpNinv);
//  P1.setBlock(1,1,Iy);
  

  TSFLinearOperator prec = P1*P2*P3;

  // TSFLinearOperator prec = P1;

  return new GenericRightPreconditioner(prec);
}

string EPWBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





