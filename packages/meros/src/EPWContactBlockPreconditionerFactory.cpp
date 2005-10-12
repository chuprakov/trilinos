/// @HEADER
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
// License along with this libxrary; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "EPWContactBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "GenericLeftPreconditioner.h"
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
#include "TSFDiagonalOperator.h"


using namespace Meros;
using namespace TSF;
using std::string;

//Need to update to reflect the Schur solver (pick a default solver)
EPWContactBlockPreconditionerFactory
::EPWContactBlockPreconditionerFactory(const TSFLinearSolver& blockSolverA,
				       const TSFLinearSolver& blockSolverC,
				       const TSFLinearSolver& blockSolverE,
				       const TSFLinearSolver& schurSolver)
  : TSFPreconditionerFactoryBase(), 
    blockSolverA_(blockSolverA),
    blockSolverC_(blockSolverC),
    blockSolverE_(blockSolverE),
    schurSolver_(schurSolver)
{}


// //Need to update to reflect the Schur solver (pick a default solver)
// EPWContactBlockPreconditionerFactory
// ::EPWContactBlockPreconditionerFactory(const TSFLinearSolver& blockSolver,
// 				const TSFLinearSolver& SchurSolver)
//  	: TSFPreconditionerFactoryBase(), 
// 	  blockSolver_(blockSolver),
// 	  SchurSolver_(SchurSolver)
// {}


TSFPreconditioner EPWContactBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  TSFError::raise("EPWContactBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFOperatorSource instead of TSFLinearOperator");
  return 0;
}

TSFPreconditioner EPWContactBlockPreconditionerFactory
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
  cerr << "problem size (domain) = " << op.domain().dim() << endl;
  cerr << "problem size (range) = " << op.range().dim() << endl;

  // Check how many blocks we have 
  int blocks = op.numBlockRows();

  TSFLinearOperator A = op.getBlock(0,0);
  TSFLinearOperator C = op.getBlock(1,1);
  TSFLinearOperator E = op.getBlock(2,2);

  //  TSFLinearOperator Kt1 = op.getBlock(0,1);
  //  TSFLinearOperator K = op.getBlock(1,0);


  TSFLinearOperator Ainv = A.inverse(blockSolverA_);
  TSFLinearOperator Cinv = C.inverse(blockSolverC_);
  TSFLinearOperator Einv = E.inverse(blockSolverE_);

  //  TSFLinearOperator Ainv = new TSFIdentityOperator(A.domain());
  //  TSFLinearOperator Cinv = new TSFIdentityOperator(C.domain());
  //  TSFLinearOperator Einv = new TSFIdentityOperator(E.domain());

  
  TSFLinearOperator Q = new TSFBlockLinearOperator(op.domain(),op.range());

  // Put the pieces together into the composed block preconditioner
  Q.setBlock(0,0,Ainv);
  Q.setBlock(1,1,Cinv);
  Q.setBlock(2,2,Einv);

  // Now make the Schur complement operator
  TSFLinearOperator Kt1 = op.getBlock(0,3);
  TSFLinearOperator Kt2 = op.getBlock(1,3);
  TSFLinearOperator Kt3 = op.getBlock(2,3);
  TSFLinearOperator K1 = op.getBlock(3,0);
  TSFLinearOperator K2 = op.getBlock(3,1);
  TSFLinearOperator K3 = op.getBlock(3,2);

  TSFLinearOperator X = K1*Ainv*Kt1 + K2*Cinv*Kt2 + K3*Einv*Kt3;
  TSFLinearOperator Xinv = X.inverse(schurSolver_);


  //  TSFLinearOperator Xinv = new TSFIdentityOperator(Kt1.domain());
  // TSFLinearOperator X = K*P*Kt;
  // TSFLinearOperator X = K*Kt;

  Q.setBlock(3,3,Xinv);

  TSFVector x = Q.domain().createMember();
  x.randomize();
  TSFVector y = Q.range().createMember();
  y.randomize();

  y = Q*x;
  cerr << "norm2(Q*x) = " << y.norm2() << endl;
  //  Q.apply(x,y);

  cerr << "Got here -- next step is return from createPreconditioner " << endl;

  return new GenericRightPreconditioner(Q);
}

string EPWContactBlockPreconditionerFactory::toString() const 
{
  return "NS block preco factory";
}





