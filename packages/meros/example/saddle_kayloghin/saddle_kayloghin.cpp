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

  cerr << "starting saddle_kayloghin" << endl;

  // Build epetra pressure and velocity maps for a simple example

  int Nvelocity = 480;
  int Npressure = 256;
  Epetra_Map *VelocityMap = new Epetra_Map(Nvelocity, 0, *comm);
  Epetra_Map *PressureMap = new Epetra_Map(Npressure, 0, *comm);
  Epetra_Map *BlockMap    = new Epetra_Map(Nvelocity + Npressure, 0, *comm);

  // Read F,B, and B^T from files and store as epetra crs matrices

  Epetra_CrsMatrix *F_crs, *Bt_crs, *B_crs;
  ReadPetraMatrix(VelocityMap, VelocityMap,  &F_crs, "../data/mac/F.serial");
  ReadPetraMatrix(VelocityMap, PressureMap, &Bt_crs, "../data/mac/Bt.serial");
  ReadPetraMatrix(PressureMap, VelocityMap,  &B_crs, "../data/mac/B.serial");

  // Build TSF vector spaces. This could probably go into a subroutine so
  // we could hide the details!!!

  Epetra_Map    *vel1 = (Epetra_Map *) &(F_crs->OperatorDomainMap());
  Epetra_Map    *vel2 = (Epetra_Map *) &(F_crs->RowMatrixColMap());
  Epetra_Import *vel3  = (Epetra_Import *) F_crs->RowMatrixImporter();
  TSFSmartPtr<Epetra_Map>    tv1 = TSFSmartPtr<Epetra_Map>(vel1,false);
  TSFSmartPtr<Epetra_Map>    tv2 = TSFSmartPtr<Epetra_Map>(vel2,false);
  TSFSmartPtr<Epetra_Import> tv3 = TSFSmartPtr<Epetra_Import>(vel3,false);
  TSFVectorSpace velocitySpace = new PetraVectorSpace( tv1, tv2, tv3);

  Epetra_Map *p1 = (Epetra_Map *) &(Bt_crs->OperatorDomainMap());
  Epetra_Map *p2 = (Epetra_Map *) &(Bt_crs->RowMatrixColMap());
  Epetra_Import *p3 = (Epetra_Import *) Bt_crs->RowMatrixImporter();
  TSFSmartPtr<Epetra_Map> tp1 = TSFSmartPtr<Epetra_Map>(p1,false);
  TSFSmartPtr<Epetra_Map> tp2 = TSFSmartPtr<Epetra_Map>(p2,false);
  TSFSmartPtr<Epetra_Import> tp3 = TSFSmartPtr<Epetra_Import>(p3,false);
  TSFVectorSpace pressureSpace = new PetraVectorSpace( tp1, tp2, tp3);


  // Convert block matrices to TSF matrices. Each individual conversion
  // could be shoved into a function to hide the details.

  PetraMatrix*  F_petra = new PetraMatrix(velocitySpace, velocitySpace);
  F_petra->setPetraMatrix(F_crs,true);   // insert epetra matrix into TSF matrix.
  TSFLinearOperator F_tsf = F_petra;     // 'true' indicates that when this TSF matrix
                                         // is deleted (automatically via smart
                                         // pointers),  the underlying epetra matrix is 
                                         // also deleted.

  PetraMatrix* Bt_petra = new PetraMatrix(velocitySpace, pressureSpace);
  Bt_petra->setPetraMatrix(Bt_crs,true); // insert epetra matrix into TSF matrix.
  TSFLinearOperator Bt_tsf = Bt_petra;   // 'true' indicates that when this TSF matrix
                                         // is deleted (automatically via smart
                                         // pointers),  the underlying epetra matrix is 
                                         // also deleted.

  PetraMatrix* B_petra = new PetraMatrix(pressureSpace, velocitySpace);
  B_petra->setPetraMatrix(B_crs,true); // insert epetra matrix into TSF matrix.
  TSFLinearOperator B_tsf = B_petra;   // 'true' indicates that when this TSF matrix
                                         // is deleted (automatically via smart
                                         // pointers),  the underlying epetra matrix is 
                                         // also deleted.

  // Set up 2x2 block TSF operator

  TSFVectorSpace  rangeBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);
  TSFVectorSpace domainBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);
  TSFLinearOperator   saddleA_tsf = new TSFBlockLinearOperator(rangeBlockSpace, domainBlockSpace);

  saddleA_tsf.setBlock(0,0,F_tsf);
  saddleA_tsf.setBlock(0,1,Bt_tsf);
  saddleA_tsf.setBlock(1,0,B_tsf);

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


 // 1) Build inv(F) so that it corresponds to using GMRES with ML.

  TSFLinearSolver FSolver;
  ML_solverData   Fsolver_data;
  bool symmetric = false;
  ML_TSF_defaults(FSolver, &Fsolver_data, symmetric, F_crs);
  FSolver.setVerbosityLevel(4);
  TSFLinearOperator F_inv = F_tsf.inverse(FSolver);

 // 2) Build a Schur complement factory for getting inv(X) approximation.

 // 2 a) Build solver for inv(Ap) 
 //      using TSF's GMRES solver
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

 TSFLinearOperator Fp_tsf = new TSFIdentityOperator(pressureSpace);
 TSFLinearOperator Ap_tsf = new TSFIdentityOperator(pressureSpace);
 TSFOperatorSource opSrc = new KayLoghinRightOperatorSource(saddleA_tsf, Fp_tsf, Ap_tsf);

 // 4 b) Create the Kay & Loghin style preconditioner as a TSFLinearOperator
 TSFPreconditioner P = pfac.createPreconditioner(opSrc);
 TSFLinearOperator saddleM_tsf = P.right();


 // Build a map. This map is intended to go between vbr-style vectors and TSF-style vectors. 
 // Since this example does not contain VBR data, the map is pointless but needed
 // for the proper creation of an epetra wrapper around a TSF matrix.
 // Specifically, map[i] != 1 indicates that the ith variable is in the
 //                           first block.
 //               map[i] == 1 indicates that the ith variable is in the
 //                           second block.

 int *imap = (int *) malloc(sizeof(int)*(Nvelocity + Npressure));
 for (int i = 0; i < Nvelocity; i++) imap[i] = 0;
 for (int i = Nvelocity; i < Nvelocity + Npressure; i++)  imap[i] = 1;

 // Build the epetra matrix wrapper around the TSF block matrix

 TSFLinearOperator2EpetraRowMatrix *saddleA_epet = 
   new TSFLinearOperator2EpetraRowMatrix(saddleA_tsf, comm, BlockMap, imap,
					 EPETRA_MATRIX);

 // Build the epetra matrix wrapper around the TSF block PRECONDITIONER

 TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet 
   = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf, comm, BlockMap,
	               (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                       (saddleA_epet))->getBlockAssignments(), EPETRA_INVERSE);
  
 // Read in 'x' and 'rhs'
 
 Epetra_Vector *x_epet   = new Epetra_Vector(*BlockMap);
 Epetra_Vector *rhs_epet = new Epetra_Vector(*BlockMap);
 ReadPetraVector(x_epet  , "../data/mac/init_guess");
 ReadPetraVector(rhs_epet, "../data/mac/rhs");
 

 // Set up the outer iteration
 
 Epetra_LinearProblem problem(saddleA_epet, x_epet, rhs_epet);
 AztecOO solver(problem);
 solver.SetPrecOperator(saddlePrec_epet);
 solver.SetAztecOption(AZ_solver,AZ_GMRESR);
 
 solver.Iterate(30, 1.0E-8);
 
 delete saddlePrec_epet;
 delete x_epet;
 delete rhs_epet;

#ifndef EPETRA_MPI
 delete comm;
#endif
 
 
 TSFMPI::finalize(); 
}



