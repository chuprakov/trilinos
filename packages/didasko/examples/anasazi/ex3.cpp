
// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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

// Compute the lowest eigenvalues and corresponding eigenvectors using
// Anasazi::BlockKrylovSchur

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_ANASAZI) && defined(HAVE_DIDASKO_TEUCHOS) && defined(HAVE_DIDASKO_TRIUTILS) 

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_MultiVector.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"

int main(int argc, char *argv[]) {
    
#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  // If we aren't using MPI, then setup a serial communicator.
  Epetra_SerialComm Comm;
#endif

  // Some typedefs for oft-used data types
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

  int ierr;

  // Get our processor ID (0 if running serial)
  int MyPID = Comm.MyPID();

  // Verbose flag: only processor 0 should print
  bool verbose = (MyPID==0);

  // Initialize a Gallery object, from which we will select a matrix
  // Select recirc_2d, of order 100
  Trilinos_Util::CrsMatrixGallery Gallery("recirc_2d", Comm);
  Gallery.Set("problem_size", 100);

  // Say hello and print some information about the gallery matrix
  if (verbose) {
    cout << "Anasazi Example: Block Kyrlov Schur" << endl;
    cout << "Problem info:" << endl;
    cout << Gallery;
    cout << endl;
  }

  // Setup some more problem/solver parameters:

  // Number of desired eigenvalues
  int nev = 4;

  // Block size
  int blocksize = 4;


  // Get a pointer to the system matrix, inside of a Teuchos::RefCountPtr
  // The Teuchos::RefCountPtr is a reference counting pointer that handles
  // garbage collection for us, so that we can perform memory allocation without
  // having to worry about freeing memory manually.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Gallery.GetMatrix() );


  // Create an Anasazi MultiVector, based on Epetra MultiVector
  const Epetra_Map * Map = &(A->RowMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = 
    Teuchos::rcp( new Epetra_MultiVector(*Map,blocksize) );

  // Fill it with random numbers
  ivec->Random();


  // Setup the eigenproblem, with the matrix A and the initial vectors ivec
  Teuchos::RefCountPtr< Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(A, ivec) );

  // The 2-D laplacian is symmetric. Specify this in the eigenproblem.
  MyProblem->SetSymmetric(false);

  // Specify the desired number of eigenvalues
  MyProblem->SetNEV( nev );

  // Signal that we are done setting up the eigenvalue problem
  ierr = MyProblem->SetProblem();

  // Check the return from SetProblem(). If this is non-zero, there was an
  // error. This probably means we did not specify enough information for
  // the eigenproblem.
  assert(ierr==0);


  // Create an output manager
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = 
    Teuchos::rcp(  new Anasazi::OutputManager<double>( MyPID ) );

  // Specify the verbosity level. Options include:
  // Anasazi::Error 
  //   This option is always set
  // Anasazi::Warning 
  //   Warnings (less severe than errors)
  // Anasazi::IterationDetails 
  //   Details at each iteration, such as the current eigenvalues
  // Anasazi::OrthoDetails 
  //   Details about orthogonality
  // Anasazi::TimingDetails
  //   A summary of the timing info for the solve() routine
  // Anasazi::FinalSummary 
  //   A final summary (we will call this manually below)
  // Anasazi::Debug 
  //   Debugging information
  MyOM->SetVerbosity( Anasazi::Warning );


  // Create the parameter list for the eigensolver
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blocksize);
  MyPL.set( "Max Blocks", 5 );
  MyPL.set( "Max Restarts", 300 );
  MyPL.set( "Tol", 1.0e-8 );

  // Create the sort manager
  // Choices are:
  // LM - target the largest magnitude  [default]
  // SM - target the smallest magnitude 
  // LR - target the largest real 
  // SR - target the smallest real 
  // LI - target the largest imaginary
  // SI - target the smallest imaginary
  std::string which("SM");
  Teuchos::RefCountPtr<Anasazi::BasicSort<double,MV,OP> > MySM =
    Teuchos::rcp( new Anasazi::BasicSort<double,MV,OP>(which) );

  // Create the Block Krylov Schur solver
  // This takes as inputs the eigenvalue problem, the sort manager
  // that we specified above, the output manager, and the solver parameters
  Anasazi::BlockKrylovSchur<double,MV,OP> 
    MyBlockKrylovSchur(MyProblem, MySM, MyOM, MyPL );

  // Solve the eigenvalue problem, and save the return code
  Anasazi::ReturnType solverreturn = MyBlockKrylovSchur.solve();

  // Print results from eigensolver
  MyBlockKrylovSchur.currentStatus();

  // Check return code of the solver: Unconverged, Failed, or OK
  switch (solverreturn) {

  // UNCONVERGED
  case Anasazi::Unconverged:
    if (verbose) 
      cout << "Anasazi::BlockKrylovSchur::solve() did not converge!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    break;

  // FAILED
  case Anasazi::Failed:
    if (verbose)
      cout << "Anasazi::BlockKrylovSchur::solve() failed!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
    break;

  // OK
  case Anasazi::Ok:
    if (verbose) 
      cout << "Anasazi::BlockKrylovSchur::solve()\n" 
           << "        converged in " << MyBlockKrylovSchur.GetNumIters() 
           << " iterations with " << MyBlockKrylovSchur.GetNumRestarts() 
           << " restarts." << endl;
    break;
  }

  // Get eigenvalues
  // Because the matrix is nonsymmetric (and the eigenvalues are potentially 
  // complex), GetEvals() returns a vector of length 2*nev:
  // evals = [ real(lambda_1)   ]
  //         [      ...
  //         [ real(lambda_nev) ]
  //         [ imag(lambda_1)   ]
  //         [      ...
  //         [ imag(labmda_nev) ]
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();

  // Get eigenvectors
  // This returns a multivector, with 2*nev vectors
  // The real part of the eigenvectors is in the first nev vectors
  // The imaginary part of the eigenvectors is in the last nev vectors
  Teuchos::RefCountPtr< Epetra_MultiVector > evecs = MyProblem->GetEvecs();
  // Extract the real and imaginary parts of the eigenvectors into two
  // separate vectors
  Teuchos::RefCountPtr< Epetra_MultiVector > evecr = 
    Teuchos::rcp( new Epetra_MultiVector(View,*evecs,0  ,nev) );
  Teuchos::RefCountPtr< Epetra_MultiVector > eveci = 
    Teuchos::rcp( new Epetra_MultiVector(View,*evecs,nev,nev) );

  // Test residuals
  // Generate (nev x nev) dense matrices for the eigenvalues
  // These matrices are automatically initialized to zero
  Teuchos::SerialDenseMatrix<int, double> Dr(nev, nev), Di(nev, nev);

  // Put the eigenvalues on the diagonals
  for (int i=0; i<nev; i++) {
    Dr(i,i) = (*evals)[i];
    Di(i,i) = (*evals)[nev+i];
  }

  // Generate a multivector for the residual
  Epetra_MultiVector R(*Map,nev);

  // Also declare a vector<double> to hold the norm
  std::vector<double> normRr(nev), normRi(nev);

  // Compute the residual, manually performing the complex arithmetic
  // R = A*evecs - evecs*D = A*(evecr+i*eveci) - (evecr+i*eveci)*(Dr+i*Di)
  //   = (A*evecr-evecr*Dr+eveci*Di) + i*(A*eveci-evecr*Di-eveci*Dr)

  // Real part: A*evecr-evecr*Dr+eveci*Di
  A->Apply( *evecr, R );
  MVT::MvTimesMatAddMv( -1.0, *evecr, Dr, 1.0, R );
  MVT::MvTimesMatAddMv( +1.0, *eveci, Di, 1.0, R );
  
  // Compute the 2-norm of the real part of the residual vectors
  MVT::MvNorm( R, &normRr );

  // Imag part: A*eveci-evecr*Di-eveci*Dr
  A->Apply( *eveci, R );
  MVT::MvTimesMatAddMv( -1.0, *evecr, Di, 1.0, R );
  MVT::MvTimesMatAddMv( -1.0, *eveci, Dr, 1.0, R );

  // Compute the 2-norm of the imag part of the residual vectors
  MVT::MvNorm( R, &normRi );

  // We will use the LAPACK wrapper in Teuchos to compute the
  // complex absolute value
  Teuchos::LAPACK<int,double> lapack;

  for (int i=0; i<nev; i++) {
    normRr[i] = lapack.LAPY2(normRr[i],normRi[i]);
  }

  // Output results to screen
  if(verbose) {
    cout << scientific << showpos << setprecision(6) << showpoint;
    cout << "******************************************************\n"
         << "           Results (outside of eigensolver)           \n" 
         << "------------------------------------------------------\n"
         << "  real(evals)\t\timag(evals)\tnormR\n"
         << "------------------------------------------------------\n";
    for(int i=0; i<nev; ++i) {
      cout << "  " << (*evals)[i] << "\t\t" << (*evals)[nev+i] 
           << "\t" << normRr[i] << endl;
    }
    cout << "******************************************************\n" << endl;
  }


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("Please configure Didasko with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-triutils");
  puts("--enable-anasazi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif

