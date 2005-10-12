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
EPWBlockPreconditionerFactory
::EPWBlockPreconditionerFactory(const TSFLinearSolver& blockSolverA,
				const TSFLinearSolver& blockSolverC,
				const TSFLinearSolver& blockSolverE)
 	: TSFPreconditionerFactoryBase(), 
	  blockSolverA_(blockSolverA),
	  blockSolverC_(blockSolverC),
	  blockSolverE_(blockSolverE)
{}

// //Need to update to reflect the Schur solver (pick a default solver)
// EPWBlockPreconditionerFactory
// ::EPWBlockPreconditionerFactory(const TSFLinearSolver& blockSolver,
// 				const TSFLinearSolver& SchurSolver)
//  	: TSFPreconditionerFactoryBase(), 
// 	  blockSolver_(blockSolver),
// 	  SchurSolver_(SchurSolver)
// {}


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

  // VEH 3/24/2005
  // original matrix is
  //      |    A    B    |
  //      |    B    C    |
  // simpler block solver:
  //      |  inv(A)   0     |
  //      |    0    inv(C)  |
  // another simple block solver:
  //      |  inv(A)   0     |
  //      |    0    inv(X)  |
  //  X = approx Schur complement

  // Check how many blocks we have 
  int blocks = op.numBlockRows();
  cerr << "num blocks for prec is " << blocks << endl;


  TSFLinearOperator A = op.getBlock(0,0);
  // TSFLinearOperator B = op.getBlock(0,1);
  TSFLinearOperator C = op.getBlock(1,1);
  // TSFLinearOperator Dinv = new TSFDiagonalOperator;
  // Epetra_CrsMatrix *C_crs = PetraMatrix::getConcrete(C);
  // cerr << "num diags = " << C_crs->NumGlobalDiagonals() << endl;

  // Epetra_Vector Cdiags(C_crs->Map());
  // C_crs->ExtractDiagonalCopy(Cdiags);
  // get the reciprocals
  //  Epetra_Vector CdiagsInv(C_crs->Map());
  //  CdiagsInv.Reciprocal(Cdiags);

  TSFLinearOperator Ainv = A.inverse(blockSolverA_);
  TSFLinearOperator Cinv = C.inverse(blockSolverC_);

  TSFLinearOperator E = op.getBlock(2,2);
  TSFLinearOperator Einv = E.inverse(blockSolverE_);
  

  // Make block operators for the three preconditioner factors
  TSFLinearOperator P = new TSFBlockLinearOperator(op.domain(), op.range());

  // Put the pieces together into the composed block preconditioner
  P.setBlock(0,0,Ainv);
  P.setBlock(1,1,Cinv);

  P.setBlock(2,2,Einv);

  TSFLinearOperator prec = P;

  //  return new GenericRightPreconditioner(prec);
  return new GenericLeftPreconditioner(prec);
}

string EPWBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





