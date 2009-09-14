/*! \@HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas\@sandia.gov)

************************************************************************
*/
/*! \@HEADER */


#include "CTrilinos_config.h"

#include "CEpetra_Map.h"
#include "CEpetra_Time.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_JadMatrix.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Flops.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Operator.h"
#include "CEpetra_DistObject.h"
#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_CompObject.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#include "mpi.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "../../../../epetra/test/epetra_test_err.h"
#include "Epetra_Version.h"
#include <vector>
#include <algorithm>
#include <string>

// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    return(1); 
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}

int checkMultiVectors( CT_Epetra_MultiVector_ID_t & X, CT_Epetra_MultiVector_ID_t & Y, string message = "", bool verbose = false) {
  int numVectors = Epetra_MultiVector_NumVectors(X);
  int length = Epetra_MultiVector_MyLength(Y);
  int badvalue = 0;
  int globalbadvalue = 0;
  for (int j=0; j<numVectors; j++) {
    CT_Epetra_Vector_ID_t vecx = Epetra_MultiVector_getVector(X, j);
    CT_Epetra_Vector_ID_t vecy = Epetra_MultiVector_getVector(Y, j);
    for (int i=0; i< length; i++)
      if (checkValues(Epetra_Vector_getElement(vecx, i), Epetra_Vector_getElement(vecy, i))==1) badvalue = 1;
    Epetra_Vector_Destroy(&vecy);
    Epetra_Vector_Destroy(&vecx);
  }
  CT_Epetra_DistObject_ID_t doX = Epetra_DistObject_Cast(Epetra_MultiVector_Abstract(X));
  CT_Epetra_Comm_ID_t Comm = Epetra_DistObject_Comm(doX);
  Epetra_Comm_MaxAll_Int(Comm, &badvalue, &globalbadvalue, 1);
  Epetra_Comm_Destroy(&Comm);
  Epetra_DistObject_Destroy(&doX);

  if (verbose) {
    if (globalbadvalue==0) cout << message << " check OK." << endl;
    else cout << "********* " << message << " check failed.********** " << endl;
  }
  return(globalbadvalue);
}

int check(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & B, bool verbose);

int powerMethodTests(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & JadA, CT_Epetra_Map_ID_t & Map, 
		     CT_Epetra_Vector_ID_t & q, CT_Epetra_Vector_ID_t & z, CT_Epetra_Vector_ID_t & resid, bool verbose);

int power_method(bool TransA, CT_Epetra_RowMatrix_ID_t& A, 
		 CT_Epetra_Vector_ID_t& q,
		 CT_Epetra_Vector_ID_t& z0, 
		 CT_Epetra_Vector_ID_t& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int main(int argc, char *argv[])
{
  int ierr = 0, i, forierr = 0;
#ifdef HAVE_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  CT_Epetra_Comm_ID_t Comm = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(Epetra_MpiComm_Create( MPI_COMM_WORLD )));

#else

  int rank = 0;
  CT_Epetra_Comm_ID_t Comm = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(Epetra_SerialComm_Create()));

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Epetra_Comm_Broadcast_Int(Comm, &verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Epetra_Comm_Barrier(Comm);

//  Epetra_Comm_SetTracebackMode(Comm, 0); // This should shut down any error traceback reporting
  int MyPID = Epetra_Comm_MyPID(Comm);
  int NumProc = Epetra_Comm_NumProc(Comm);

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) 
		verbose = false;

  int NumMyEquations = 10000;
  int NumGlobalEquations = (NumMyEquations * NumProc) + EPETRA_MIN(NumProc,3);
  if(MyPID < 3) 
    NumMyEquations++;

  // Construct a Map that puts approximately the same Number of equations on each processor

  CT_Epetra_Map_ID_t Map = Epetra_Map_Create_Linear(NumGlobalEquations, NumMyEquations, 0, Comm);
  CT_Epetra_BlockMap_ID_t bMap = Epetra_BlockMap_Cast(Epetra_Map_Abstract(Map));
  
  // Get update list and number of local equations from newly created Map
  vector<int> MyGlobalElements(Epetra_BlockMap_NumMyElements(bMap));
  Epetra_BlockMap_MyGlobalElements_Fill(bMap, &MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  vector<int> NumNz(NumMyEquations);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for(i = 0; i < NumMyEquations; i++)
    if((MyGlobalElements[i] == 0) || (MyGlobalElements[i] == NumGlobalEquations - 1))
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_Matrix

  CT_Epetra_CrsMatrix_ID_t A = Epetra_CrsMatrix_Create_VarPerRow(
      CT_Epetra_DataAccess_E_Copy, Map, &NumNz[0], FALSE);
  EPETRA_TEST_ERR(Epetra_CrsMatrix_IndicesAreGlobal(A),ierr);
  EPETRA_TEST_ERR(Epetra_CrsMatrix_IndicesAreLocal(A),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  vector<double> Values(2);
  Values[0] = -1.0; 
	Values[1] = -1.0;
	vector<int> Indices(2);
  double two = 2.0;
  int NumEntries;

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) {
    if(MyGlobalElements[i] == 0) {
			Indices[0] = 1;
			NumEntries = 1;
		}
    else if (MyGlobalElements[i] == NumGlobalEquations-1) {
			Indices[0] = NumGlobalEquations-2;
			NumEntries = 1;
		}
    else {
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			NumEntries = 2;
		}
		forierr += !(Epetra_CrsMatrix_InsertGlobalValues(A, MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		forierr += !(Epetra_CrsMatrix_InsertGlobalValues(A, MyGlobalElements[i], 1, &two, &MyGlobalElements[i])>0); // Put in the diagonal entry
  }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  Epetra_CrsMatrix_FillComplete(A, TRUE);
  Epetra_CrsMatrix_OptimizeStorage(A);

  CT_Epetra_RowMatrix_ID_t rA = Epetra_RowMatrix_Cast(Epetra_CrsMatrix_Abstract(A));
  CT_Epetra_JadMatrix_ID_t JadA = Epetra_JadMatrix_Create(rA);
  CT_Epetra_RowMatrix_ID_t rJadA = Epetra_RowMatrix_Cast(Epetra_JadMatrix_Abstract(JadA));
  CT_Epetra_JadMatrix_ID_t JadA1 = Epetra_JadMatrix_Create(rA);
  CT_Epetra_RowMatrix_ID_t rJadA1 = Epetra_RowMatrix_Cast(Epetra_JadMatrix_Abstract(JadA1));
  CT_Epetra_JadMatrix_ID_t JadA2 = Epetra_JadMatrix_Create(rA);
  CT_Epetra_RowMatrix_ID_t rJadA2 = Epetra_RowMatrix_Cast(Epetra_JadMatrix_Abstract(JadA2));

  // Create vectors for Power method

  CT_Epetra_Vector_ID_t q = Epetra_Vector_Create(bMap, TRUE);
  CT_Epetra_Vector_ID_t z = Epetra_Vector_Create(bMap, TRUE);
  CT_Epetra_MultiVector_ID_t mz = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(z));
  Epetra_MultiVector_Random(mz);
  CT_Epetra_Vector_ID_t resid = Epetra_Vector_Create(bMap, TRUE);

  CT_Epetra_Flops_ID_t flopcounter = Epetra_Flops_Create();

  CT_Epetra_CompObject_ID_t coA = Epetra_CompObject_Cast(Epetra_CrsMatrix_Abstract(A));
  Epetra_CompObject_SetFlopCounter(coA, flopcounter);

  CT_Epetra_CompObject_ID_t coq = Epetra_CompObject_Cast(Epetra_Vector_Abstract(q));
  Epetra_CompObject_SetFlopCounter_Matching(coq, coA);

  CT_Epetra_CompObject_ID_t coz = Epetra_CompObject_Cast(Epetra_Vector_Abstract(z));
  Epetra_CompObject_SetFlopCounter_Matching(coz, coA);

  CT_Epetra_CompObject_ID_t coresid = Epetra_CompObject_Cast(Epetra_Vector_Abstract(resid));
  Epetra_CompObject_SetFlopCounter_Matching(coresid, coA);

  CT_Epetra_CompObject_ID_t coJadA = Epetra_CompObject_Cast(Epetra_JadMatrix_Abstract(JadA));
  Epetra_CompObject_SetFlopCounter_Matching(coJadA, coA);

  CT_Epetra_CompObject_ID_t coJadA1 = Epetra_CompObject_Cast(Epetra_JadMatrix_Abstract(JadA1));
  Epetra_CompObject_SetFlopCounter_Matching(coJadA1, coA);

  CT_Epetra_CompObject_ID_t coJadA2 = Epetra_CompObject_Cast(Epetra_JadMatrix_Abstract(JadA2));
  Epetra_CompObject_SetFlopCounter_Matching(coJadA2, coA);
  

  if (verbose) cout << "=======================================" << endl
		    << "Testing Jad using CrsMatrix as input..." << endl
		    << "=======================================" << endl;

  Epetra_CompObject_ResetFlops(coA);
  powerMethodTests(rA, rJadA, Map, q, z, resid, verbose);

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (Epetra_CrsMatrix_MyGlobalRow(A, 0)) {
    int numvals = Epetra_CrsMatrix_NumGlobalEntries(A, 0);
    vector<double> Rowvals(numvals);
    vector<int> Rowinds(numvals);
    Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices(A, 0, numvals, &numvals, &Rowvals[0], &Rowinds[0]); // Get A[0,0]

    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
    
    Epetra_CrsMatrix_ReplaceGlobalValues(A, 0, numvals, &Rowvals[0], &Rowinds[0]);
  }
  Epetra_JadMatrix_UpdateValues(JadA, rA, FALSE);
  Epetra_CompObject_ResetFlops(coA);
  powerMethodTests(rA, rJadA, Map, q, z, resid, verbose);

  if (verbose) cout << "================================================================" << endl
		          << "Testing Jad using Jad matrix as input matrix for construction..." << endl
		          << "================================================================" << endl;
  Epetra_CompObject_ResetFlops(coJadA1);
  powerMethodTests(rJadA1, rJadA2, Map, q, z, resid, verbose);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}

int powerMethodTests(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & JadA, CT_Epetra_Map_ID_t & Map, 
		     CT_Epetra_Vector_ID_t & q, CT_Epetra_Vector_ID_t & z, CT_Epetra_Vector_ID_t & resid, bool verbose) {

  // variable needed for iteration
  double lambda = 0.0;
  // int niters = 10000;
  int niters = 300;
  double tolerance = 1.0e-2;
  int ierr = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Iterate

  CT_Epetra_BlockMap_ID_t bMap = Epetra_BlockMap_Cast(Epetra_Map_Abstract(Map));
  CT_Epetra_Comm_ID_t Comm = Epetra_BlockMap_Comm(bMap);
  Epetra_BlockMap_Destroy(&bMap);
  CT_Epetra_Time_ID_t timer = Epetra_Time_Create(Comm);
  Epetra_Comm_Destroy(&Comm);
	
  CT_Epetra_CompObject_ID_t coq = Epetra_CompObject_Cast(Epetra_Vector_Abstract(q));

  double startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  double elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  double total_flops = Epetra_CompObject_Flops(coq);
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  double lambdaref = lambda;
  double flopsref = total_flops;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for reference first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(false, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(coq);
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref," No-transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref," No-transpose Power Method flop count", verbose),ierr);

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;
  // Iterate
  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(coq);
  MFLOPs = total_flops/elapsed_time/1000000.0;
  lambdaref = lambda;
  flopsref = total_flops;

  if (verbose) 
	 cout << "\n\nTotal MFLOPs for reference transpose solve = " << MFLOPs << endl
		 <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(true, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(coq);
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate transpose solve = " << MFLOPs << endl
		  <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref,"Transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref,"Transpose Power Method flop count", verbose),ierr);

  EPETRA_TEST_ERR(check(A, JadA, verbose),ierr);

  Epetra_CompObject_Destroy(&coq);

  return(0);
}
int power_method(bool TransA, CT_Epetra_RowMatrix_ID_t& A, CT_Epetra_Vector_ID_t& q, CT_Epetra_Vector_ID_t& z0, 
		 CT_Epetra_Vector_ID_t& resid, double* lambda, int niters, double tolerance, bool verbose) 
{  
	
  // Fill z with random Numbers
  CT_Epetra_Vector_ID_t z = Epetra_Vector_Duplicate(z0);
	
  CT_Epetra_MultiVector_ID_t mz = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(z));
  CT_Epetra_MultiVector_ID_t mq = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(q));
  CT_Epetra_MultiVector_ID_t mresid = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(resid));

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;
	
  for(int iter = 0; iter < niters; iter++) {
		Epetra_MultiVector_Norm2(mz, &normz); // Compute 2-norm of z
		Epetra_MultiVector_Scale(mq, 1.0/normz, mz);
		Epetra_RowMatrix_Multiply(A, TransA, mq, mz); // Compute z = A*q // SEGFAULT HAPPENS HERE
		Epetra_MultiVector_Dot(mq, mz, lambda); // Approximate maximum eigenvaluE
		if(iter%100==0 || iter+1==niters) {
			Epetra_MultiVector_Update_WithAB(mresid, 1.0, mz, -(*lambda), mq, 0.0); // Compute A*q - lambda*q
			Epetra_MultiVector_Norm2(mresid, &residual);
			if(verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
											 << "  Residual of A*q - lambda*q = " << residual << endl;
		}
		if(residual < tolerance) {
			ierr = 0;
			break;
		}
	}

  Epetra_MultiVector_Destroy(&mresid);
  Epetra_MultiVector_Destroy(&mq);
  Epetra_MultiVector_Destroy(&mz);

  Epetra_Vector_Destroy(&z);

  return(ierr);
}

int check(CT_Epetra_RowMatrix_ID_t& A, CT_Epetra_RowMatrix_ID_t & B, bool verbose)  {  

  int ierr = 0;

  CT_Epetra_Operator_ID_t oA = Epetra_Operator_Cast(Epetra_RowMatrix_Abstract(A));
  CT_Epetra_Operator_ID_t oB = Epetra_Operator_Cast(Epetra_RowMatrix_Abstract(B));

  CT_Epetra_Comm_ID_t CommA = Epetra_Operator_Comm(oA);
  CT_Epetra_Comm_ID_t CommB = Epetra_Operator_Comm(oB);
  EPETRA_TEST_ERR(!Epetra_Comm_NumProc(CommA)==Epetra_Comm_NumProc(CommB),ierr);
  EPETRA_TEST_ERR(!Epetra_Comm_MyPID(CommA)==Epetra_Comm_MyPID(CommB),ierr);
  Epetra_Comm_Destroy(&CommB);
  Epetra_Comm_Destroy(&CommA);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_Filled(A)==Epetra_RowMatrix_Filled(B),ierr);
  EPETRA_TEST_ERR(!Epetra_Operator_HasNormInf(oA)==Epetra_Operator_HasNormInf(oB),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_LowerTriangular(A)==Epetra_RowMatrix_LowerTriangular(B),ierr);

  CT_Epetra_SrcDistObject_ID_t sdoA = Epetra_SrcDistObject_Cast(Epetra_RowMatrix_Abstract(A));
  CT_Epetra_SrcDistObject_ID_t sdoB = Epetra_SrcDistObject_Cast(Epetra_RowMatrix_Abstract(B));
  CT_Epetra_BlockMap_ID_t mapA = Epetra_SrcDistObject_Map(sdoA);
  CT_Epetra_BlockMap_ID_t mapB = Epetra_SrcDistObject_Map(sdoB);
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(mapA, mapB),ierr);
  Epetra_BlockMap_Destroy(&mapB);
  Epetra_BlockMap_Destroy(&mapA);
  Epetra_SrcDistObject_Destroy(&sdoB);
  Epetra_SrcDistObject_Destroy(&sdoA);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_MaxNumEntries(A)==Epetra_RowMatrix_MaxNumEntries(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalCols(A)==Epetra_RowMatrix_NumGlobalCols(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalDiagonals(A)==Epetra_RowMatrix_NumGlobalDiagonals(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalNonzeros(A)==Epetra_RowMatrix_NumGlobalNonzeros(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalRows(A)==Epetra_RowMatrix_NumGlobalRows(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyCols(A)==Epetra_RowMatrix_NumMyCols(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyDiagonals(A)==Epetra_RowMatrix_NumMyDiagonals(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyNonzeros(A)==Epetra_RowMatrix_NumMyNonzeros(B),ierr);
  for (int i=0; i<Epetra_RowMatrix_NumMyRows(A); i++) {
    int nA, nB;
    Epetra_RowMatrix_NumMyRowEntries(A,i,&nA);
    Epetra_RowMatrix_NumMyRowEntries(B,i,&nB);
    EPETRA_TEST_ERR(!nA==nB,ierr);
  }
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyRows(A)==Epetra_RowMatrix_NumMyRows(B),ierr);

  CT_Epetra_Map_ID_t odmA = Epetra_Operator_OperatorDomainMap(oA);
  CT_Epetra_BlockMap_ID_t bodmA = Epetra_BlockMap_Cast(Epetra_Map_Abstract(odmA));
  Epetra_Map_Destroy(&odmA);
  CT_Epetra_Map_ID_t odmB = Epetra_Operator_OperatorDomainMap(oB);
  CT_Epetra_BlockMap_ID_t bodmB = Epetra_BlockMap_Cast(Epetra_Map_Abstract(odmB));
  Epetra_Map_Destroy(&odmB);
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(bodmA, bodmB),ierr);

  CT_Epetra_Map_ID_t ormA = Epetra_Operator_OperatorRangeMap(oA);
  CT_Epetra_BlockMap_ID_t bormA = Epetra_BlockMap_Cast(Epetra_Map_Abstract(ormA));
  Epetra_Map_Destroy(&ormA);
  CT_Epetra_Map_ID_t ormB = Epetra_Operator_OperatorRangeMap(oB);
  CT_Epetra_BlockMap_ID_t bormB = Epetra_BlockMap_Cast(Epetra_Map_Abstract(ormB));
  Epetra_Map_Destroy(&ormB);
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(bormA, bormB),ierr);

  CT_Epetra_Map_ID_t rcmA = Epetra_RowMatrix_RowMatrixColMap(A);
  CT_Epetra_BlockMap_ID_t brcmA = Epetra_BlockMap_Cast(Epetra_Map_Abstract(rcmA));
  Epetra_Map_Destroy(&rcmA);
  CT_Epetra_Map_ID_t rcmB = Epetra_RowMatrix_RowMatrixColMap(B);
  CT_Epetra_BlockMap_ID_t brcmB = Epetra_BlockMap_Cast(Epetra_Map_Abstract(rcmB));
  Epetra_Map_Destroy(&rcmB);
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(brcmA, brcmB),ierr);

  CT_Epetra_Map_ID_t rrmA = Epetra_RowMatrix_RowMatrixRowMap(A);
  CT_Epetra_BlockMap_ID_t brrmA = Epetra_BlockMap_Cast(Epetra_Map_Abstract(rrmA));
  Epetra_Map_Destroy(&rrmA);
  CT_Epetra_Map_ID_t rrmB = Epetra_RowMatrix_RowMatrixRowMap(B);
  CT_Epetra_BlockMap_ID_t brrmB = Epetra_BlockMap_Cast(Epetra_Map_Abstract(rrmB));
  Epetra_Map_Destroy(&rrmB);
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(brrmA, brrmB),ierr);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_UpperTriangular(A)==Epetra_RowMatrix_UpperTriangular(B),ierr);
  EPETRA_TEST_ERR(!Epetra_Operator_UseTranspose(oA)==Epetra_Operator_UseTranspose(oB),ierr);

  int NumVectors = 5;
  { // No transpose case
    CT_Epetra_MultiVector_ID_t X = Epetra_MultiVector_Create(bodmA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA1 = Epetra_MultiVector_Create(bormA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA2 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB1 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB2 = Epetra_MultiVector_Duplicate(YA1);
    Epetra_MultiVector_Random(X);
    
    bool transA = false;
    Epetra_Operator_SetUseTranspose(oA, transA);
    Epetra_Operator_SetUseTranspose(oB, transA);
    Epetra_Operator_Apply(oA,X,YA1);
    Epetra_RowMatrix_Multiply(A, transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2,"A Multiply and A Apply", verbose),ierr);
    Epetra_Operator_Apply(oB,X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1,"A Multiply and B Multiply", verbose),ierr);
    Epetra_RowMatrix_Multiply(B, transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2,"A Multiply and B Apply", verbose), ierr);

    Epetra_MultiVector_Destroy(&X);
    Epetra_MultiVector_Destroy(&YB2);
    Epetra_MultiVector_Destroy(&YB1);
    Epetra_MultiVector_Destroy(&YA2);
    Epetra_MultiVector_Destroy(&YA1);
    
  }
  {// transpose case
    CT_Epetra_MultiVector_ID_t X = Epetra_MultiVector_Create(bormA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA1 = Epetra_MultiVector_Create(bodmA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA2 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB1 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB2 = Epetra_MultiVector_Duplicate(YA1);
    Epetra_MultiVector_Random(X);
    
    bool transA = true;
    Epetra_Operator_SetUseTranspose(oA, transA);
    Epetra_Operator_SetUseTranspose(oB, transA);
    Epetra_Operator_Apply(oA,X,YA1);
    Epetra_RowMatrix_Multiply(A, transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2, "A Multiply and A Apply (transpose)", verbose),ierr);
    Epetra_Operator_Apply(oB,X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1, "A Multiply and B Multiply (transpose)", verbose),ierr);
    Epetra_RowMatrix_Multiply(B, transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2, "A Multiply and B Apply (transpose)", verbose),ierr);
    
    Epetra_MultiVector_Destroy(&X);
    Epetra_MultiVector_Destroy(&YB2);
    Epetra_MultiVector_Destroy(&YB1);
    Epetra_MultiVector_Destroy(&YA2);
    Epetra_MultiVector_Destroy(&YA1);
    
  }

  Epetra_Operator_Destroy(&oB);
  Epetra_Operator_Destroy(&oA);

  CT_Epetra_Vector_ID_t diagA = Epetra_Vector_Create(brrmA, TRUE);
  CT_Epetra_MultiVector_ID_t mdiagA = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(diagA));
  EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractDiagonalCopy(A, diagA),ierr);
  CT_Epetra_Vector_ID_t diagB = Epetra_Vector_Create(brrmB, TRUE);
  CT_Epetra_MultiVector_ID_t mdiagB = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(diagB));
  EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractDiagonalCopy(B, diagB),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(mdiagA,mdiagB, "ExtractDiagonalCopy", verbose),ierr);

  CT_Epetra_Vector_ID_t rowA = Epetra_Vector_Create(brrmA, TRUE);
  CT_Epetra_MultiVector_ID_t mrowA = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(rowA));
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvRowSums(A, rowA),ierr);
  CT_Epetra_Vector_ID_t rowB = Epetra_Vector_Create(brrmB, TRUE);
  CT_Epetra_MultiVector_ID_t mrowB = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(rowB));
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvRowSums(B, rowB),ierr)
  EPETRA_TEST_ERR(checkMultiVectors(mrowA,mrowB, "InvRowSums", verbose),ierr);

  CT_Epetra_Vector_ID_t colA = Epetra_Vector_Create(brcmA, TRUE);
  CT_Epetra_MultiVector_ID_t mcolA = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(colA));
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvColSums(A, colA),ierr);
  CT_Epetra_Vector_ID_t colB = Epetra_Vector_Create(brcmB, TRUE);
  CT_Epetra_MultiVector_ID_t mcolB = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(colB));
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvColSums(B, colB),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(mcolA,mcolB, "InvColSums", verbose),ierr);

  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormInf(A), Epetra_RowMatrix_NormInf(B), "NormInf before scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormOne(A), Epetra_RowMatrix_NormOne(B), "NormOne before scaling", verbose),ierr);

  EPETRA_TEST_ERR(Epetra_RowMatrix_RightScale(A, colA),ierr);  
  EPETRA_TEST_ERR(Epetra_RowMatrix_RightScale(B, colB),ierr);


  EPETRA_TEST_ERR(Epetra_RowMatrix_LeftScale(A, rowA),ierr);
  EPETRA_TEST_ERR(Epetra_RowMatrix_LeftScale(B, rowB),ierr);

  Epetra_MultiVector_Destroy(&mcolB);
  Epetra_MultiVector_Destroy(&mcolA);
  Epetra_MultiVector_Destroy(&mrowB);
  Epetra_MultiVector_Destroy(&mrowA);
  Epetra_MultiVector_Destroy(&mdiagB);
  Epetra_MultiVector_Destroy(&mdiagA);

  Epetra_Vector_Destroy(&colB);
  Epetra_Vector_Destroy(&colA);
  Epetra_Vector_Destroy(&rowB);
  Epetra_Vector_Destroy(&rowA);
  Epetra_Vector_Destroy(&diagB);
  Epetra_Vector_Destroy(&diagA);

  Epetra_BlockMap_Destroy(&brrmA);
  Epetra_BlockMap_Destroy(&brrmB);
  Epetra_BlockMap_Destroy(&brcmA);
  Epetra_BlockMap_Destroy(&brcmB);
  Epetra_BlockMap_Destroy(&bormA);
  Epetra_BlockMap_Destroy(&bormB);
  Epetra_BlockMap_Destroy(&bodmA);
  Epetra_BlockMap_Destroy(&bodmB);

  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormInf(A), Epetra_RowMatrix_NormInf(B), "NormInf after scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormOne(A), Epetra_RowMatrix_NormOne(B), "NormOne after scaling", verbose),ierr);

  vector<double> valuesA(Epetra_RowMatrix_MaxNumEntries(A));
  vector<int> indicesA(Epetra_RowMatrix_MaxNumEntries(A));  
  vector<double> valuesB(Epetra_RowMatrix_MaxNumEntries(B));
  vector<int> indicesB(Epetra_RowMatrix_MaxNumEntries(B));
  return(0);
  for (int i=0; i<Epetra_RowMatrix_NumMyRows(A); i++) {
    int nA, nB;
    EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractMyRowCopy(A, i, Epetra_RowMatrix_MaxNumEntries(A), &nA, &valuesA[0], &indicesA[0]),ierr); 
    EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractMyRowCopy(B, i, Epetra_RowMatrix_MaxNumEntries(B), &nB, &valuesB[0], &indicesB[0]),ierr);
    EPETRA_TEST_ERR(!nA==nB,ierr);
    for (int j=0; j<nA; j++) {
      double curVal = valuesA[j];
      int curIndex = indicesA[j];
      bool notfound = true;
      int jj = 0;
      while (notfound && jj< nB) {
	if (!checkValues(curVal, valuesB[jj])) notfound = false;
	jj++;
      }
      EPETRA_TEST_ERR(notfound, ierr);
      vector<int>::iterator p = find(indicesB.begin(),indicesB.end(),curIndex);  // find curIndex in indicesB
      EPETRA_TEST_ERR(p==indicesB.end(), ierr);
    }

  }
  if (verbose) cout << "RowMatrix Methods check OK" << endl;

  return (ierr);
}
