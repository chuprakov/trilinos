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

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Vector.h"
#include "Vbr2Petra.h"
#include "ReadPetra.h"
#include "TSFLinearOperator2EpetraRowMatrix.h"
#include "ml_include.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"
#include "TSFIdentityOperator.h"
#include "Aztec2TSF.h"

int main(int argc, char *argv[])
{
  TSFMPI::init(&argc, (void***)&argv); // (cast to void for tsf)


  // Read Matrix from file into Aztec VBR matrix
  
  int Nnz = 4228, Nblks = 256, blk_size = 3;
  AZ_MATRIX *Amat;

  Epetra_Comm *comm;
#ifdef EPETRA_MPI
  Epetra_MpiComm mpicomm(MPI_COMM_WORLD);
  comm = (Epetra_Comm *) (&mpicomm);
#else
  comm = (Epetra_Comm *) new Epetra_SerialComm();
#endif

  // Convert Aztec matrix to a set of 4 block epetra matrices

  Epetra_Map        *VbrMap;
  Epetra_Map **subblock_maps;
  subblock_maps = (Epetra_Map **) malloc(sizeof(Epetra_Map *)*2);
  subblock_maps[0] = new Epetra_Map(Nblks*(blk_size-1), 0, *comm);
  subblock_maps[1] = new Epetra_Map(Nblks, 0, *comm);
  ReadAztecVbr(blk_size, Nnz, Nblks, &Amat, "../data/tmac-vbr/K_reordered");

  Epetra_RowMatrix *saddleA_epet = Aztec2TSF(Amat,comm,(Epetra_BlockMap *&) VbrMap,subblock_maps);

  // Pull out stuff from the matrix 

  TSFLinearOperator saddleA_tsf = (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>(saddleA_epet))->getTSF();
  TSFLinearOperator  F_tsf = saddleA_tsf.getBlock(0,0);
  Epetra_CrsMatrix  *F_crs = PetraMatrix::getConcrete(F_tsf);
  Epetra_Map        *F_map = (Epetra_Map *) &(F_crs->OperatorDomainMap());

  TSFLinearOperator  B_tsf = saddleA_tsf.getBlock(1,0);
  TSFLinearOperator  Bt_tsf = saddleA_tsf.getBlock(0,1);
  TSFLinearOperator BBt_tsf;
  TSF_MatrixMult(B_tsf, Bt_tsf, BBt_tsf);

  Epetra_CrsMatrix  *BBt_crs = PetraMatrix::getConcrete(BBt_tsf);
  Epetra_Map        *BBt_map = (Epetra_Map *) &(BBt_crs->OperatorDomainMap());



 // build a simple preconditioner corresponding to 
 //
 //                    |   inv(F)     0    |
 //                    |     0        I    |
 //
 // 1) Build inv(F) so that it corresponds to using GMRES with ML.


  TSFLinearSolver FSolver,      BBtSolver;
  ML_solverData   Fsolver_data, BBtSolver_data;

  bool symmetric = false;

  ML_TSF_defaults(FSolver, &Fsolver_data, symmetric, F_crs);
  FSolver.setVerbosityLevel(4);
  TSFLinearOperator F_inv = F_tsf.inverse(FSolver);

  symmetric = true;

  ML_TSF_defaults(BBtSolver, &BBtSolver_data, symmetric, BBt_crs);
  BBtSolver.setVerbosityLevel(4);
  TSFLinearOperator BBt_inv = BBt_tsf.inverse(BBtSolver);

   // 2) Build 2x2 block preconditioner by seting inv(F) to the (0,0)
   //    block and the identity to the (1,1) block. Then wrap the 
   //    preconditioner within an Epetra_Row matrix.

   TSFLinearOperator saddleM_tsf =  new TSFBlockLinearOperator(saddleA_tsf.range(), saddleA_tsf.domain());
   saddleM_tsf.setBlock(0,0,F_inv);
   TSFLinearOperator ident = new TSFIdentityOperator(saddleA_tsf.getBlock(1,0).range());
   //      saddleM_tsf.setBlock(1,1,ident);
       saddleM_tsf.setBlock(1,1,BBt_inv);

   TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf,
                                             comm, (Epetra_Map *) &(saddleA_epet->OperatorDomainMap()),
(dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>(saddleA_epet))->getBlockAssignments(),
					     EPETRA_INVERSE);

   // Read in 'x' and 'rhs'

   Epetra_Vector *x_epet, *rhs_epet;
   x_epet   = new Epetra_Vector(*VbrMap);
   rhs_epet = new Epetra_Vector(*VbrMap);
   ReadPetraVector(x_epet  , "../data/tmac-vbr/guess_reordered");
   ReadPetraVector(rhs_epet, "../data/tmac-vbr/rhs_reordered");

   // Set up the outer iteration

   Epetra_LinearProblem problem(saddleA_epet, x_epet, rhs_epet);
   AztecOO solver(problem);
   solver.SetPrecOperator(saddlePrec_epet);
   solver.SetAztecOption(AZ_solver,AZ_GMRESR);
   solver.Iterate(200, 1.0E-8);




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



