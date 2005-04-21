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
#include "TSFPreconditioner.h"
#include "TSFMatrixOperator.h"
#include "Epetra_Vector.h"
#include "TSFHashtable.h"

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
#include "KayLoghinSchurFactory.h"
#include "GMRESSolver.h"
#include "Aztec2TSF.h"

#include "ml_include.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"

#include "TSFIdentityOperator.h"

using namespace TSF;
using namespace Meros;

int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  Epetra_Comm *comm;
#ifdef EPETRA_MPI
  Epetra_MpiComm mpicomm(MPI_COMM_WORLD);
  comm = (Epetra_Comm *) (&mpicomm);
#else
  comm = (Epetra_Comm *) new Epetra_SerialComm();
#endif

  // Read Matrix from file into Aztec VBR matrix
  
  cerr << "starting saddle_kayloghin" << endl;
  int Nnz = 4228, Nblks = 256, blk_size = 3;
  AZ_MATRIX *Amat;
  ReadAztecVbr(blk_size, Nnz, Nblks, &Amat, "../data/mac-vbr/K_reordered");

  // Convert Aztec matrix to a set of 4 block epetra matrices

  Epetra_Map         *VbrMap;
  Epetra_Map         **subblock_maps;
  subblock_maps    = (Epetra_Map **) malloc(sizeof(Epetra_Map *)*2);
  subblock_maps[0] = new Epetra_Map(Nblks*(blk_size-1), 0, *comm);
  subblock_maps[1] = new Epetra_Map(Nblks, 0, *comm);

  Epetra_RowMatrix *saddleA_epet=Aztec2TSF(Amat,comm,(Epetra_BlockMap *&) VbrMap,subblock_maps);

  // Pull out stuff from the matrix 

  TSFLinearOperator saddleA_tsf = (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                                  (saddleA_epet))->getTSF();
  TSFLinearOperator F_tsf = saddleA_tsf.getBlock(0,0);
  TSFLinearOperator B_tsf = saddleA_tsf.getBlock(1,0);
  Epetra_CrsMatrix *F_crs = PetraMatrix::getConcrete(F_tsf);
  Epetra_Map       *F_map = (Epetra_Map *) &(F_crs->OperatorDomainMap());
  TSFVectorSpace   pressureSpace = B_tsf.range();

 // Build a Kay & Loghin style block preconditioner with meros
 // 
 //      |  inv(F)   0   | |   I    -Bt  | |   I      0     |
 //      |    0      I   | |   0     I   | |   0   -inv(X)  |
 //

 // We'll do this in 4 steps.
 // 1) Build a solver for inv(F)
 // 2) Build a SchurFactory that can make an inv(X) approximation
 // 3) Build a block preconditioner factory with the F solver and Schur factory
 // 4) Make the preconditioner and get a TSFLinearOperator representing the prec. 


 // 1) Build solver for inv(F) so that it corresponds to using GMRES with ML.
  TSFHashtable<int, int> azOptionsF;
  TSFHashtable<int, double> azParamsF;
  azOptionsF.put(AZ_solver, AZ_gmres);
  azOptionsF.put(AZ_ml, 1);
  azOptionsF.put(AZ_ml_levels, 10);
  azOptionsF.put(AZ_precond, AZ_dom_decomp);
  azOptionsF.put(AZ_subdomain_solve, AZ_ilu);
  azParamsF.put(AZ_tol, 1e-6);
  azOptionsF.put(AZ_max_iter, 200);
  azOptionsF.put(AZ_recursive_iterate, 1);
  TSFLinearSolver FSolver = new AZTECSolver(azOptionsF,azParamsF);
  FSolver.setVerbosityLevel(4); 

  // int N_levels = 10;
  // ML_Set_PrintLevel(10);
  // ML_Create(&ml_handle, N_levels);
  // EpetraMatrix2MLMatrix(ml_handle, 0, F_crs);
  // ML_Aggregate *agg_object;
  // ML_Aggregate_Create(&agg_object);
  // ML_Aggregate_Set_MaxCoarseSize(agg_object,30);
  // ML_Set_Symmetrize(ml_handle, ML_TRUE);
  // ML_Aggregate_Set_DampingFactor(agg_object,0.25);
  // N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
  //                               ML_INCREASING, agg_object);
  // ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS,
  //                              ML_BOTH, 2, .2);
  // ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);

  // Epetra_ML_Operator  *MLop = 
  //                 new Epetra_ML_Operator(ml_handle,*comm,*F_map,*F_map);
  // MLop->SetOwnership(true);
  // ML_Aggregate_Destroy(&agg_object);
  // TSFSmartPtr<Epetra_Operator> Smart_MLprec = 
  //                    TSFSmartPtr<Epetra_Operator>(MLop, true);

  //  printf("the way ML is passed to an aztec solver has changed\n");
  // printf("so this code no longer works. Someone will have to fix\n");
  // printf("this\n");
  // exit(1);
  // TSFLinearSolver FSolver;//= new AZTECSolver(azOptions, azParams, Smart_MLprec);
 FSolver.setVerbosityLevel(4);
 
  TSFLinearOperator F_inv = F_tsf.inverse(FSolver);

 
 // 2) Build a Schur complement factory for getting inv(X) approximation.

 // 2 a) Build solver for inv(Ap) so that it corresponds to using CG with ML.
 //      Using same MLop as built for F solver.
 // Can't use aztec for Ap solve yet since my Ap isn't an epetra matrix.
 // azOptions.put(AZ_solver, AZ_cg);
 // azOptions.put(AZ_conv, AZ_r0);
 // azParams.put(AZ_tol, 1e-8);
 // azOptions.put(AZ_max_iter, 250);
 // azOptions.put(AZ_precond, AZ_none);
 // TSFLinearSolver ApSolver = new AZTECSolver(azOptions, azParams, &MLop);
 // TSFLinearSolver ApSolver = new AZTECSolver(azOptions, azParams);
 // TSF's unrestarted GMRES
 TSFLinearSolver ApSolver = new GMRESSolver(1e-08, 250, 250);
 ApSolver.setVerbosityLevel(1);

 // 2 b) Build a Schur complement factory of type KayLoghinSchurFactory.
 SchurFactory sfac = new KayLoghinSchurFactory(ApSolver);
 
 // 3) Build a preconditioner factory with the F solver and the Schur factory.
 TSFPreconditionerFactory pfac = new NSBlockPreconditionerFactory(FSolver, sfac);

 // 4) Give the preconditioner factory data (the linear operators) 
 //    and build the preconditioner.

 // 4 a) Build a Kay and Loghin type operator source.
 //      The operator source contains all of the linear ops needed for the preconditioner.
 //      For a Kay and Loghin preconditioner we need saddleA, Fp, and Ap.
 //      Setting Fp and Ap to identity operators for now. Later will read these in.

 // Make an epetra identity matrix so I can try Aztec/ML in our preconditioner Ap solve
 //  PetraMatrix* Ap_petra = new PetraMatrix(pressureSpace, pressureSpace);
 // putting in C since it's the right size
 //  for (int index = 0; index < pressureSpace.dim(); index++)
 //    Ap_petra->setElement(index, index, diag);
//  double val = 1.0;
//  for (int index = 0; index < pressureSpace.dim(); index++) 
//    C->InsertGlobalValues(index,1,&val,&index);
//  C->FillComplete();
//  cout << "got here" << endl;
//  Ap_petra->setPetraMatrix(C);    // insert the epetra matrix into our TSF 
//  TSFLinearOperator Ap_tsf = Ap_petra; 
//  TSFReal diag = 1.0;

//  Epetra_CrsMatrix *Fp_crs;
//  ReadPetraMatrix(&p1, &Fp_crs, "FpMatrixFile");
//  Epetra_CrsMatrix *Ap_crs;
//  ReadPetraMatrix(&p1, &Ap_crs, "ApMatrixFile");

//  PetraMatrix* Fp_petra = new PetraMatrix(pressureSpace, pressureSpace);
//  Fp_petra->setPetraMatrix(Fp_crs);    // insert the epetra matrix into our TSF 
//  TSFLinearOperator Fp_tsf = Fp_petra; // matrix and make it a TSF linear op
//  PetraMatrix* Ap_petra = new PetraMatrix(pressureSpace, pressureSpace);
//  Ap_petra->setPetraMatrix(Ap_crs);    // insert the epetra matrix into our TSF 
//  TSFLinearOperator Ap_tsf = Ap_petra; // matrix and make it a TSF linear op

 TSFLinearOperator Fp_tsf = new TSFIdentityOperator(pressureSpace);
 TSFLinearOperator Ap_tsf = new TSFIdentityOperator(pressureSpace);
 TSFOperatorSource opSrc = new KayLoghinRightOperatorSource(saddleA_tsf, Fp_tsf, Ap_tsf);

 // 4 b) Create the Kay & Loghin style preconditioner as a TSFLinearOperator
 TSFPreconditioner P = pfac.createPreconditioner(opSrc);
 TSFLinearOperator saddleM_tsf = P.right();

//  TSFVector xtmp = domainBlockSpace.createMember();
//  TSFVector ytmp = rangeBlockSpace.createMember();

//  saddleM_tsf.apply(xtmp, ytmp);

//  exit(0);



 // Convert the TSFLinearOperator to an EpetraRowMatrix
 TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet 
   = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf, comm, VbrMap, 
	               (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                       (saddleA_epet))->getBlockAssignments(), EPETRA_INVERSE);
  
 // Read in 'x' and 'rhs'
 
 Epetra_Vector *x_epet, *rhs_epet;
 x_epet   = new Epetra_Vector(*VbrMap);
 rhs_epet = new Epetra_Vector(*VbrMap);
 ReadPetraVector(x_epet  , "../data/mac-vbr/guess_reordered");
 ReadPetraVector(rhs_epet, "../data/mac-vbr/rhs_reordered");
 
 // Set up the outer iteration
 
 
 Epetra_LinearProblem problem(saddleA_epet, x_epet, rhs_epet);
 AztecOO solver(problem);
 solver.SetPrecOperator(saddlePrec_epet);
 solver.SetAztecOption(AZ_solver,AZ_GMRESR);
 
 solver.Iterate(3, 1.0E-8);
 
 
 delete saddlePrec_epet;
 delete saddleA_epet;
 delete x_epet;
 delete rhs_epet;
 delete subblock_maps[0];
 delete subblock_maps[1];

#ifndef EPETRA_MPI
 delete comm;
#endif

 delete VbrMap;
 
 
 free(Amat->data_org);
 free(Amat->bpntr);
 free(Amat->bindx);
 free(Amat->val);
 free(Amat->cpntr);
 free(Amat->rpntr);
 free(Amat->indx);
 free(subblock_maps);
 AZ_matrix_destroy(&Amat);
 
 
 TSFMPI::finalize(); 
}



