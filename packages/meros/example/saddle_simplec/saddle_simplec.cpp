#include "simplec.h"

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

  cerr << "starting pressure projection" << endl;

  // Build epetra pressure and velocity maps for a simple example
  int Nvelocity = 578;
  int Npressure = 256;
  Epetra_Map *VelocityMap = new Epetra_Map(Nvelocity, 0, *comm);
  Epetra_Map *PressureMap = new Epetra_Map(Npressure, 0, *comm);
  Epetra_Map *BlockMap    = new Epetra_Map(Nvelocity + Npressure, 0, *comm);

  // Read F,B, and B^T from files and store as epetra crs matrices
  Epetra_CrsMatrix *F_crs, *Bt_crs, *B_crs, *C_crs;
  Epetra_Vector *x_epet   = new Epetra_Vector(*BlockMap),*rhs_epet = new Epetra_Vector(*BlockMap);
 //  ReadPetraMatrix(VelocityMap, VelocityMap,  &F_crs, "../data/q1/Aq1");
//   ReadPetraMatrix(VelocityMap, PressureMap, &Bt_crs, "../data/q1/Btq1");
//   ReadPetraMatrix(PressureMap, VelocityMap,  &B_crs, "../data/q1/Bq1");
//   //  ReadPetraMatrix(PressureMap, PressureMap,  &C_crs, "../data/q1p0/C");
//   ReadPetraVector(x_epet  , "../data/q1/newrhs");
//   ReadPetraVector(rhs_epet, "../data/q1/rhsq1");

  ReadPetraMatrix(VelocityMap, VelocityMap,  &F_crs, "../data/q1p0/F");
  ReadPetraMatrix(VelocityMap, PressureMap, &Bt_crs, "../data/q1p0/Bt");
  ReadPetraMatrix(PressureMap, VelocityMap,  &B_crs, "../data/q1p0/B");
  ReadPetraMatrix(PressureMap, PressureMap,  &C_crs, "../data/q1p0/C");
  ReadPetraVector(x_epet  , "../data/q1p0/initial");
  ReadPetraVector(rhs_epet, "../data/q1p0/rhs");
  
  // Build TSF vector spaces. This could probably go into a subroutine so
  // we could hide the details!!!
  TSFVectorSpace velocitySpace = EpetraCRS2TSFVspace(F_crs);
  TSFVectorSpace pressureSpace = EpetraCRS2TSFVspace(Bt_crs);

  // Convert block matrices to TSF matrices. Each individual conversion
  // could be shoved into a function to hide the details.
  TSFLinearOperator F_tsf  = EpetraCRS2TSF(velocitySpace,velocitySpace,F_crs);
  TSFLinearOperator Bt_tsf = EpetraCRS2TSF(pressureSpace,velocitySpace,Bt_crs);
  TSFLinearOperator B_tsf  = EpetraCRS2TSF(velocitySpace,pressureSpace,B_crs);
  TSFLinearOperator C_tsf  = EpetraCRS2TSF(pressureSpace, pressureSpace, C_crs);  
  //  C_tsf.describe();
  // Set up 2x2 block TSF operator
  TSFVectorSpace  rangeBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);
  TSFVectorSpace domainBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);
  TSFLinearOperator   saddleA_tsf = new TSFBlockLinearOperator(rangeBlockSpace, domainBlockSpace);

  saddleA_tsf.setBlock(0,0,F_tsf);
  saddleA_tsf.setBlock(0,1,Bt_tsf);
  saddleA_tsf.setBlock(1,0,B_tsf);
  saddleA_tsf.setBlock(1,1,C_tsf);
  saddleA_tsf.describe();
 // Build a pressure projection style block preconditioner with meros
 // 
 //      |  I  -inv(diag(F))Bt | |   F       0     |-1
 //      |  0          I       | |   B      inv(X) |
 // where inv(X) = -Binv(diag(F))Bt.
 // We'll do this in 4 steps:
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

 // 2) Build a Schur complement factory for getting inv(X) approximation.
 // 2 a) Build solver for the Schur complement using TSF's GMRES solver
  /*  double tol = 1e-8;
   int maxiter = 250;
   TSFHashtable<int, int>azOptions;
   TSFHashtable<int, double>azParams;
   azOptions.put(AZ_solver, AZ_cg);
   azOptions.put(AZ_precond, AZ_none);
   azParams.put(AZ_tol, tol);
   azOptions.put(AZ_max_iter, maxiter);
   TSFLinearSolver SchurSolver = new AZTECSolver(azOptions, azParams);
  */
   TSFLinearSolver SchurSolver = new GMRESSolver(1e-8,250,250);
   SchurSolver.setVerbosityLevel(4);

 // 2 b) Build a Schur complement factory of SIMPLE type.
   SchurFactory sfac = new SimpleSchurFactory(SchurSolver);
 
 // 3) Build a preconditioner factory with the F solver and the Schur factory.
 TSFPreconditionerFactory pfac = new SimpleCBlockPreconditionerFactory(FSolver, sfac, SchurSolver);

 // 4) Give the preconditioner factory data (the linear operators) 
 //    and build the preconditioner.
 // 4 a) Build a SIMPLE type operator source.
 //      The operator source contains all of the linear ops needed for the preconditioner.
 //      For a SIMPLE preconditioner we need saddleA.
  TSFOperatorSource opSrc = new SimpleCOperatorSource(saddleA_tsf);

 // 4 b) Create the SIMPLE style preconditioner as a TSFLinearOperator
 TSFPreconditioner P = pfac.createPreconditioner(opSrc);
 TSFLinearOperator saddleM_tsf = P.right();

 // This map is intended to go between vbr-style vectors and TSF-style vectors. 
 // The map is needed for the proper creation of an epetra wrapper around a TSF matrix.
  int *imap;
  imap = Vbr2TSF(Nvelocity, Npressure);

 // Build the epetra matrix wrapper around the TSF block matrix
 TSFLinearOperator2EpetraRowMatrix *saddleA_epet = 
   new TSFLinearOperator2EpetraRowMatrix(saddleA_tsf, comm, BlockMap, imap,
					 EPETRA_MATRIX);

 // Build the epetra matrix wrapper around the TSF block PRECONDITIONER
 TSFLinearOperator2EpetraRowMatrix *saddlePrec_epet 
   = new TSFLinearOperator2EpetraRowMatrix(saddleM_tsf, comm, BlockMap,
	               (dynamic_cast<TSFLinearOperator2EpetraRowMatrix*>
                       (saddleA_epet))->getBlockAssignments(), EPETRA_INVERSE);
  
 // Set up the outer iteration 
 Epetra_LinearProblem problem(saddleA_epet, x_epet, rhs_epet);
 AztecOO solver(problem);
 solver.SetPrecOperator(saddlePrec_epet);
 solver.SetAztecOption(AZ_solver,AZ_GMRESR);
 solver.Iterate(290, 1.0E-6);
 
 delete saddlePrec_epet;
 delete x_epet;
 delete rhs_epet;
#ifndef EPETRA_MPI
 delete comm;
#endif
 TSFMPI::finalize(); 
}
