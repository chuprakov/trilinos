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
  ReadAztecVbr(blk_size, Nnz, Nblks, &Amat, "../data/mac-vbr/K_reordered");

  Epetra_RowMatrix *saddleA_epet = Aztec2TSF(Amat,comm,(Epetra_BlockMap *&) VbrMap,subblock_maps);

  // Pull out stuff from the matrix 

  TSFLinearOperator saddleA_tsf = (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>(saddleA_epet))->getTSF();
  TSFLinearOperator  F_tsf = saddleA_tsf.getBlock(0,0);
  Epetra_CrsMatrix  *F_crs = PetraMatrix::getConcrete(F_tsf);
  Epetra_Map        *F_map = (Epetra_Map *) &(F_crs->OperatorDomainMap());


 // build a simple preconditioner corresponding to 
 //
 //                    |   inv(F)     0    |
 //                    |     0        I    |
 //
 // 1) Build inv(F) so that it corresponds to using GMRES with ML.

   ML *ml_handle;
   int N_levels = 10;
   ML_Set_PrintLevel(10);
   ML_Create(&ml_handle, N_levels);
   EpetraMatrix2MLMatrix(ml_handle, 0, F_crs);
   ML_Aggregate *agg_object;
   ML_Aggregate_Create(&agg_object);
   ML_Aggregate_Set_MaxCoarseSize(agg_object,30);
   ML_Set_Symmetrize(ml_handle, ML_TRUE);
   ML_Aggregate_Set_DampingFactor(agg_object,0.25);
   N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                  ML_INCREASING, agg_object);
   ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 2, .2);
   ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);

   Epetra_ML_Operator  *MLop = new Epetra_ML_Operator(ml_handle,*comm,*F_map,*F_map);
   MLop->SetOwnership(true);
   ML_Aggregate_Destroy(&agg_object);


   TSFHashtable<int, int> azOptions;
   TSFHashtable<int, double> azParams;
   azOptions.put(AZ_solver, AZ_gmres);
   azOptions.put(AZ_kspace, 50);
   azOptions.put(AZ_conv, AZ_r0);
   azParams.put(AZ_tol, 1e-8);
   azOptions.put(AZ_max_iter, 50);
   TSFSmartPtr<Epetra_Operator> Smart_MLprec = TSFSmartPtr<Epetra_Operator>(MLop, true);

   TSFLinearSolver FSolver = new AZTECSolver(azOptions, azParams, Smart_MLprec);
   FSolver.setVerbosityLevel(4);
   TSFLinearOperator F_inv = F_tsf.inverse(FSolver);

   // 2) Build 2x2 block preconditioner by seting inv(F) to the (0,0)
   //    block and the identity to the (1,1) block. Then wrap the 
   //    preconditioner within an Epetra_Row matrix.

   TSFLinearOperator saddleM_tsf =  new TSFBlockLinearOperator(saddleA_tsf.range(), saddleA_tsf.domain());
   saddleM_tsf.setBlock(0,0,F_inv);
   TSFLinearOperator ident = new TSFIdentityOperator(saddleA_tsf.getBlock(1,0).range());
   saddleM_tsf.setBlock(1,1,ident);

   TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf,
                                             comm, (Epetra_Map *) &(saddleA_epet->OperatorDomainMap()),
(dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>(saddleA_epet))->getBlockAssignments(),
					     EPETRA_INVERSE);

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
   solver.Iterate(10, 1.0E-8);




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



