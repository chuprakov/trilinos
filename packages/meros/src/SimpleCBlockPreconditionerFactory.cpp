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

#include "SimpleCBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"
#include "RightBlockNSOperatorSource.h"
#include "SimpleCOperatorSource.h"
#include "BlockForwardsolver.h"
#include "TSFZeroOperator.h"
#include "Aztec2TSF.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra2TSFutils.h"
#include "Epetra_CrsMatrix.h"

using namespace Meros;
using namespace TSF;
using std::string;
//Need to update to reflect the Schur solver (pick a default solver)
SimpleCBlockPreconditionerFactory
::SimpleCBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac)
 	: TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac)
{}

SimpleCBlockPreconditionerFactory
::SimpleCBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac, const TSFLinearSolver& Schursolver)
  : TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac), Schursolver_(Schursolver)
{}

TSFPreconditioner SimpleCBlockPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& op) const
{
  TSFError::raise("SimpleBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFLinearOperator instead of TSFOperatorSource");
  return 0;
}

TSFPreconditioner SimpleCBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  // SIMPLEC style block preconditioner
  // 
  //      | I   -inv(diag(F))Bt | |   F        0               |-1
  //      | 0        I          | |   B     -B inv(diag(F)) Bt |
  //

  // Get the saddle point operator S
  TSFLinearOperator S = opSource.getOp();

   const SimpleCOperatorSource* diagSrcPtr = dynamic_cast<const SimpleCOperatorSource*>(opSource.ptr());

  // Get the blocks from S
  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Finv = diagSrcPtr->getDinv(); 
  
  TSFLinearOperator Bt = S.getBlock(0,1);
  TSFLinearOperator B = S.getBlock(1,0);
  TSFLinearOperator C = S.getBlock(1,1);
  TSFLinearOperator X, T, W;
  TSFLinearOperator FinvBt = Finv*Bt;
  
  //  TSF_MatrixMult(Finv,Bt,FinvBt);
  FinvBt.describe();
  double scal = 1.0;
  if (C.isZeroOperator()) X = B*FinvBt; //TSF_MatrixMult(B,FinvBt,X);  //DO WE WANT A MINUS C?
  else
  { 
  cerr << "Using C stabilizer" << endl;
  //TSF_MatrixMult(B, FinvBt,T);
  //TSF_MatrixAdd(C,T,scal,X);
  //X = C + B*FinvBt; 
  Epetra_CrsMatrix *Finv_crs = PetraMatrix::getConcrete(Finv);
  Epetra_CrsMatrix *Bt_crs = PetraMatrix::getConcrete(Bt);
  Epetra_CrsMatrix *B_crs = PetraMatrix::getConcrete(B);
  Epetra_CrsMatrix *C_crs = PetraMatrix::getConcrete(C);

  Epetra_CrsMatrix *BFinv_crs = new Epetra_CrsMatrix(*B_crs);
  Epetra_CrsMatrix *BFinvBt_crs = new Epetra_CrsMatrix(*C_crs);

  EpetraExt::MatrixMatrix::Multiply(*B_crs, false, *Finv_crs, false, *BFinv_crs);
  EpetraExt::MatrixMatrix::Multiply(*BFinv_crs, false, *Bt_crs, false, *BFinvBt_crs);

  EpetraExt::MatrixMatrix::Add(*C_crs, false, 1.0, *BFinvBt_crs, 1.0);

  PetraMatrix* X_petra = new PetraMatrix(C.domain(), C.range());
  X_petra->setPetraMatrix(BFinvBt_crs,true);
  X = X_petra;

  }

  // Make identity matrices on the right spaces
  TSFLinearOperator Iv = new TSFIdentityOperator(F.domain());
  TSFLinearOperator Ip = new TSFIdentityOperator(Bt.domain());

  // Make block operators for the three preconditioner factors
  TSFLinearOperator upper = new TSFBlockLinearOperator(S.domain(), S.range());
  TSFLinearOperator lower = new TSFBlockLinearOperator(S.domain(), S.range());

  // Put the pieces together into the composed block preconditioner
   lower.setBlock(0,0,F);
   lower.setBlock(1,0,B);
   lower.setBlock(1,1,X);

   upper.setBlock(0,0,Iv);
   upper.setBlock(0,1,FinvBt);
   upper.setBlock(1,1,Ip);	

   //Set up the TSF lower triangular solver with the F and Schur complement solvers
   TSFLinearSolver sol = new BlockForwardsolver(Fsolver_,Schursolver_); 
   TSFLinearOperator Linv = lower.inverse(sol);
   TSFLinearOperator prec = upper*Linv;

  return new GenericRightPreconditioner(prec);
}

string SimpleCBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





