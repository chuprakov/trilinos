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

#include <stdio.h>
#include <fstream>

#include "TSFMPI.h"
#include "TSFOut.h"
#include "TSFVectorType.h"
#include "TSFProductSpace.h"
#include "TSFBlockLinearOperator.h"
#include "PetraVectorType.h"
#include "PetraVectorSpace.h"
#include "PetraVector.h"
#include "PetraMatrix.h"
#include "GMRESSolver.h"
#include "AZTECSolver.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "ILUKRightPreconditionerFactory.h"
#include "TSFPreconditioner.h"
#include "TSFMatrixOperator.h"
#include "Epetra_Vector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Vbr2Petra.h"
#include "ReadPetra.h"
#include "TSFLinearOperator2EpetraRowMatrix.h"
#include "NSBlockPreconditionerFactory.h"
#include "SchurFactory.h"
#include "SchurFactoryBase.h"
#include "BJBSchurFactory.h"
#include "BJBRightOperatorSource.h"
#include "GMRESSolver.h"
#include "Aztec2TSF.h"

#include "ml_include.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"

#include "TSFIdentityOperator.h"

using namespace TSF;
using namespace SPP;

int main(int argc, void** argv)
{
  
  // Set up communication stuff for parallel

  TSFMPI::init(&argc, &argv);

  Epetra_Comm *comm;
#ifdef EPETRA_MPI
  Epetra_MpiComm mpicomm(MPI_COMM_WORLD);
  comm = (Epetra_Comm *) (&mpicomm);
#else
  comm = (Epetra_Comm *) new Epetra_SerialComm();
#endif
  
  cerr << "starting saddle_bjb" << endl;

  // Using salsa vbr data since it has a nonzero C block 

  // for data/mac-vbr
  //  int Nnz = 4228, Nblks = 256, blk_size = 3;
  // sizes for data/salsa
  int Nnz = 327519, Nblks = 4131, blk_size = 3;
  AZ_MATRIX *Amat;
  ReadAztecVbr(blk_size, Nnz, Nblks, &Amat, "../data/salsa/K_reordered");

  // Convert Aztec matrix to a set of 4 block epetra matrices

  Epetra_Map         *VbrMap;
  Epetra_Map         **subblock_maps;
  subblock_maps    = (Epetra_Map **) malloc(sizeof(Epetra_Map *)*2);
  subblock_maps[0] = new Epetra_Map(Nblks*(blk_size-1), 0, *comm);
  subblock_maps[1] = new Epetra_Map(Nblks, 0, *comm);

  Epetra_RowMatrix *saddleA_epet=Aztec2TSF(Amat,comm,(Epetra_BlockMap *&) VbrMap,subblock_maps);

  // Pull out stuff from the matrix 
  //      |    F    Bt   | 
  //      |    B    C    |

  TSFLinearOperator saddleA_tsf = (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                                   (saddleA_epet))->getTSF();
  TSFLinearOperator F_tsf = saddleA_tsf.getBlock(0,0);
  TSFLinearOperator Bt_tsf = saddleA_tsf.getBlock(0,1);
  TSFLinearOperator B_tsf = saddleA_tsf.getBlock(1,0);
  TSFLinearOperator C_tsf = saddleA_tsf.getBlock(1,1);

//   char Cfile[100];
//   sprintf(Cfile, "C.dat");
//   ofstream ofs(Cfile);
//   ofs.precision(15);
//   if (ofs.bad())
//     TSFError::raise("failed to open file");

  // const TSFSmartPtr<const TSFMatrixOperator> C_mat = C_tsf.getMatrix();

  //   C_mat.print(cout);

//  cerr << C_tsf << endl;


  Epetra_CrsMatrix *F_crs = PetraMatrix::getConcrete(F_tsf);
  Epetra_Map       *F_map = (Epetra_Map *) &(F_crs->OperatorDomainMap());
  TSFVectorSpace   pressureSpace = B_tsf.range();

  
  
  // Build a block preconditioner using X = inv(Ap)(B*diag(F)*Bt)inv(Ap)
  // In this case we will be using Ap = C
  // 
  //      |  inv(F)   0   | |   I    -Bt  | |   I      0     |
  //      |    0      I   | |   0     I   | |   0   -inv(X)  |
  //
  
  // We'll do this in 4 steps.
  // 1) Build a solver for inv(F)
  // 2) Build a SchurFactory that can make the inv(X) approximation
  // 3) Build a block preconditioner factory with the F solver and Schur factory
  // 4) Make the preconditioner and get a TSFLinearOperator representing the prec. 
  
  
  // 1) Build solver for inv(F) so that it corresponds to using GMRES with ML.

  TSFLinearSolver FSolver;
  ML_solverData   Fsolver_data;
  
  bool symmetric = false;

  ML_TSF_defaults(FSolver, &Fsolver_data, symmetric, F_crs);
  FSolver.setVerbosityLevel(4);
  
  TSFLinearOperator F_inv = F_tsf.inverse(FSolver);
  

  // 2) Build a Schur complement factory for getting inv(X) approximation.
  //   Using a X = inv(Ap)(B * Dinv * Bt)inv(Ap), where Ap = C;

  // solver for Xinv
  int fill = 0;
  int overlap = 0; 
  double tol = 1.0e-8;
  int maxiter = 100;
  int kspace = 100;
  // unpreconditioned GMRES
  TSFLinearSolver ApSolver = new GMRESSolver(tol, maxiter, kspace);
  ApSolver.setVerbosityLevel(4);
  

  // 2 b) Build a Schur complement factory of type BJBSchurFactory.
  SchurFactory sfac = new BJBSchurFactory(ApSolver);
 
  // 3) Build a preconditioner factory with the F solver and the Schur factory.
  TSFPreconditionerFactory pfac = new NSBlockPreconditionerFactory(FSolver, sfac);

  // 4) Give the preconditioner factory data (the linear operators) 
  //    and build the preconditioner.

  // 4 a) Build a BJBSchur type operator source.
  //      The operator source contains all of the linear ops needed for the preconditioner.
  //      For a BJBSchur preconditioner we just need saddleA.
  //      The operator source can build Dinv from A (or we can pass it in).
  
  
  TSFOperatorSource opSrc = new BJBRightOperatorSource(saddleA_tsf);


  // 4 b) Create the BJBSchur style preconditioner as a TSFLinearOperator
  TSFPreconditioner P = pfac.createPreconditioner(opSrc);
  TSFLinearOperator saddleM_tsf = P.right();
  
  // Convert the TSFLinearOperator to an EpetraRowMatrix
  TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet 
    = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf, comm, VbrMap, 
                                            (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                                             (saddleA_epet))->getBlockAssignments(), 
                                            EPETRA_INVERSE);
  
  // Read in 'x' and 'rhs'
  
  Epetra_Vector *x_epet, *rhs_epet;
  x_epet   = new Epetra_Vector(*VbrMap);
  rhs_epet = new Epetra_Vector(*VbrMap);
  ReadPetraVector(x_epet  , "../data/salsa/guess_reordered");

  // this is the wrong rhs -- there isn't one in salsa data directory
  ReadPetraVector(rhs_epet, "../data/salsa/guess_reordered");
  
  // Set up the outer iteration
 
  
  Epetra_LinearProblem problem(saddleA_epet, x_epet, rhs_epet);
  AztecOO solver(problem);
  solver.SetPrecOperator(saddlePrec_epet);
  solver.SetAztecOption(AZ_solver,AZ_GMRESR);
  
  solver.Iterate(3, 1.0E-8);
 
 

  delete saddlePrec_epet;
  delete x_epet;
  delete rhs_epet;
  
#ifndef EPETRA_MPI
  delete comm;
#endif
  
  TSFMPI::finalize(); 
}



