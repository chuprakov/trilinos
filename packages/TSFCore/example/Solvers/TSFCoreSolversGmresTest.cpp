#include "TSFCoreSolversGmresSolver.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
//
//
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	//
	// Alias MemMngPack
	namespace mmp = MemMngPack;
	//
	int i, j;
	int n_nonzeros, N_update;
	int *bindx=0, *update=0, *col_inds=0;
	double *val=0, *row_vals=0;
	double *xguess=0, *b=0, *xexact=0, *xsolve=0;

#ifdef EPETRA_MPI	
	// Initialize MPI	
	MPI_Init(&argc,&argv); 	
	Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
	Epetra_SerialComm Comm;	
#endif
	
	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	
	bool verbose = (MyPID==0);
	//
    	if(argc < 2 && verbose) {
     	cerr << "Usage: " << argv[0] 
	 << " HB_filename " << endl
	 << "where:" << endl
	 << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
	 << endl;
    	return(1);
	}
	//
	//**********************************************************************
	//******************Set up the problem to be solved*********************
	//**********************************************************************
    	//
    	int NumGlobalElements;  // total # of rows in matrix
	//
	// *****Read in matrix from HB file****** 
	//
	Trilinos_Util_read_hb(argv[1], MyPID, &NumGlobalElements, &n_nonzeros,
			      &val,  &bindx, &xguess, &b, &xexact);
	//
	// *****Distribute data among processors*****
	//
	Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update,
					 &update, &val, &bindx, &xguess, &b, &xexact);
	//
	// *****Construct the matrix*****
	//
	int NumMyElements = N_update; // # local rows of matrix on processor
	//
    	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	//
	int * NumNz = new int[NumMyElements];
	for (i=0; i<NumMyElements; i++) {
		NumNz[i] = bindx[i+1] - bindx[i] + 1;
	}
	//
	Epetra_Map Map(NumGlobalElements, NumMyElements, update, 0, Comm);
	//
	// Create a Epetra_Matrix
	//
	//Epetra_CrsMatrix A(Copy, Map, NumNz);
	mmp::ref_count_ptr<Epetra_CrsMatrix> A =
		mmp::rcp( new Epetra_CrsMatrix(Copy, Map, NumNz) );
	//
	// Create some Epetra_Vectors
	//
	mmp::ref_count_ptr<Epetra_Vector> xx = mmp::rcp( new Epetra_Vector(Copy, Map, xexact) );
	mmp::ref_count_ptr<Epetra_Vector> bb = mmp::rcp( new Epetra_Vector(Copy, Map, b) );
	// Solution vector, initialize to zero.
	mmp::ref_count_ptr<Epetra_Vector> x = mmp::rcp( new Epetra_Vector(Map) ); 
	x.get()->PutScalar( 0.0 ); 
	// 
	// Add rows one-at-a-time
	//
	int NumEntries;
	for (i=0; i<NumMyElements; i++) {
		row_vals = val + bindx[i];
		col_inds = bindx + bindx[i];
		NumEntries = bindx[i+1] - bindx[i];
		assert(A.get()->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
		assert(A.get()->InsertGlobalValues(update[i], 1, val+i, update+i)==0);
	}
	//
	// Finish up
	//
	assert(A.get()->TransformToLocal()==0);
	//
	A.get()->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
	// Create the TSFCore Linear Operator and Vectors
	//
	TSFCore::EpetraLinearOp ELOp( A ); 
	TSFCore::EpetraVector RHS( bb );
	TSFCore::EpetraVector Soln( x );
	//
    	// ********Other information used by the GMRES solver***********
	//
    	int maxits = NumGlobalElements-1; // maximum number of iterations to run
    	double tol = 1.0e-10;  // relative residual tolerance
	//
	//*******************************************************************
	// *************Start the GMRES iteration*************************
	//*******************************************************************
	//
	TSFCore::Solvers::GMRESSolver<double> MySolver( ELOp, RHS, Soln );
	//MySolver.solve();
	//
	// Compute actual residual norm.
	//
	

	if (verbose) {
	  cout << "Iteration : "<< MySolver.currIteration() << endl;
	  cout << "Final Computed GMRES Residual Norm" << 
	    MySolver.currEstRelResidualNorm() << endl;
	}

// Release all objects  

  delete [] NumNz;
  delete [] bindx;
	
  return 0;
  //
} // end TSFCoreSolversGmresTest.cpp
