
// @HEADER
// ***********************************************************************
// 
//            Amesos: An Interface to Direct Solvers
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

#include "Amesos_ConfigDefs.h"
// This example needs triutils to generate the linear system.
#ifdef HAVE_AMESOS_TRIUTILS
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util.h"

// ==================== //
// M A I N  D R I V E R //
// ==================== //
//
// This example will:
// 1.- Read an H/B matrix from file;
// 2.- redistribute the linear system matrix to the
//     available processes
// 3.- set up LHS/RHS if not present
// 4.- create an Amesos_BaseSolver object
// 5.- solve the linear problem.
//
// This example can be run with any number of processors.
//
// Author: Marzio Sala, SNL 2914
// Last modified: Apr-05

int main(int argc, char *argv[]) 
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  string matrix_file = "662_bus_out.rsa"; // file containing the HB matrix.
  if (Comm.MyPID() == 0)
    cout << "Reading matrix `" << matrix_file << "'";
  
  // ================= //
  // reading HB matrix //
  // ================= //
  
  // HB files are for serial matrices. Hence, only
  // process 0 reads this matrix (and if present
  // solution and RHS). Then, this matrix will be redistributed
  // using epetra capabilities.
  // All variables that begin with "read" refer to the
  // HB matrix read by process 0.
  Epetra_Map* readMap;
  Epetra_CrsMatrix* readA; 
  Epetra_Vector* readx; 
  Epetra_Vector* readb;
  Epetra_Vector* readxexact;
  
  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra((char*)matrix_file.c_str(), Comm, readMap, 
                              readA, readx, readb, readxexact);

  // Create uniform distributed map.
  // Note that linear map are used for simplicity only!
  // Amesos (through Epetra) can support *any* map.
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);
  
  // Create the distributed matrix, based on Map.
  Epetra_CrsMatrix A(Copy, map, 0);

  const Epetra_Map &OriginalMap = readA->RowMatrixRowMap() ; 
  assert (OriginalMap.SameAs(*readMap)); 
  Epetra_Export exporter(OriginalMap, map);
  Epetra_Export exporter2(OriginalMap, map);
  Epetra_Export MatrixExporter(OriginalMap, map);

  Epetra_Vector x(map);          // distributed solution
  Epetra_Vector b(map);          // distributed rhs
  Epetra_Vector xexact(map);     // distributed exact solution

  // Exports from data defined on processor 0 to distributed.
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);
  A.Export(*readA, exporter, Add);
  A.FillComplete();
    
  // Creates an epetra linear problem, contaning matrix
  // A, solution x and rhs b.
  Epetra_LinearProblem Problem(&A,&x,&b);
  
  // ======================================================= //
  // B E G I N N I N G   O F   T H E   A M E S O S   P A R T //
  // ======================================================= //

  // For comments on the commands in this section, please
  // see file example_AmesosFactory.cpp.
  
  string SolverType = "Klu";
  Amesos_BaseSolver* Solver = 0;
  Amesos Factory;
  
  Solver = Factory.Create(SolverType,Problem);

  // Factory.Create() returns 0 if the requested solver
  // is not available
  if (Solver == 0) {
    cerr << "Selected solver is not available" << endl;
    // return ok not to break test harness even if
    // the solver is not available
    return(EXIT_SUCCESS);
  }

  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  // =========================================== //
  // E N D   O F   T H E   A M E S O S   P A R T //
  // =========================================== //

  // Computes ||Ax - b|| //

  double residual;

  Epetra_Vector Ax(b.Map());
  A.Multiply(false, x, Ax);
  Ax.Update(1.0, b, -1.0);
  Ax.Norm2(&residual);

  if (!Comm.MyPID()) 
    cout << "After AMESOS solution, ||b-Ax||_2 = " << residual << endl;

  // delete Solver. Do this before calling MPI_Finalize() because
  // MPI calls can occur.
  delete Solver;

  if (residual > 1e-5)
    return(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

} // end of main()

#else

// Triutils is not available. Sorry, we have to give up.

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

  puts("Please configure Amesos with:");
  puts("--enable-triutils");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}

#endif
