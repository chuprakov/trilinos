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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_MPI
#include "mpi.h"
#include "CEpetra_MpiComm.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"
#include "CAmesos.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Object.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"

/*! @file example_AmesosFactory_Tridiag,c
 * This is an example of how to use the CTrilinos interface to Amesos.
 * This example follows from the Amesos example by the same name.
 */

/*
 * ====================
 * M A I N  D R I V E R
 * ====================
 *
 * This example will:
 * 1.- Create an tridiagonal matrix;
 * 2.- Call SymbolicFactorization();
 * 3.- Change the numerical values of the matrix;
 * 4.- Call NumericFactorization();
 * 5.- Set the entries of the RHS;
 * 6.- Call Solve().
 *
 * This example is intended to show the required data
 * for each phase. Phase (2) requires the matrix structure only. 
 * Phase (4) requires the matrix structure (supposed unchanged 
 * from phase (2)) and the matrix data. Phase (6) requires the 
 * RHS and solution vector.
 *
 * This example can be run with any number of processors.
 *
 * Author: Marzio Sala, SNL 9214
 * Last modified: Apr-05.
 */

int main(int argc, char *argv[]) 
{
#ifdef HAVE_CTRILINOS_AMESOS

  int NumGlobalElements, NumMyElements, NumEntries, Indices[3], i, check;
  double Values[3], residual;
  double sfact_time, nfact_time, solve_time, mtx_conv_time, mtx_redist_time, vec_redist_time;
  int *MyGlobalElements = NULL;

  CT_Epetra_Comm_ID_t Comm;
  CT_Epetra_Map_ID_t Map;
  CT_Epetra_BlockMap_ID_t bMap;
  CT_Epetra_CrsMatrix_ID_t A;
  CT_Epetra_LinearProblem_ID_t Problem;
  CT_Epetra_RowMatrix_ID_t rA;
  CT_Amesos_ID_t Factory;
  CT_Amesos_BaseSolver_ID_t Solver;
  CT_Epetra_Vector_ID_t b, x, Ax;
  CT_Epetra_MultiVector_ID_t mb, mx, mAx;
  CT_Teuchos_ParameterList_ID_t TimingsList;
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Comm = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(Epetra_MpiComm_Create(MPI_COMM_WORLD)));
#else
  Comm = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(Epetra_SerialComm_Create()));
#endif

  NumGlobalElements = 100; /* global dimension of the problem. */

  /*
   * =======================================================
   * B E G I N N I N G   O F   M A T R I X   C R E A T I O N
   * =======================================================
   */

  /*
   * Construct a Map that puts approximatively the same number of 
   * equations on each processor. `0' is the index base (that is,
   * numbering starts from 0.
   */
  Map = Epetra_Map_Create(NumGlobalElements, 0, Comm);
  bMap = Epetra_BlockMap_Cast(Epetra_Map_Abstract(Map));

  /* Create an empty EpetraCrsMatrix */
  A = Epetra_CrsMatrix_Create(CT_Epetra_DataAccess_E_Copy, Map, 0, FALSE);

  /* Create the structure of the matrix (tridiagonal) */
  NumMyElements = Epetra_BlockMap_NumMyElements(bMap);

  /* Add  rows one-at-a-time
   * Need some vectors to help
   */

  /* Right now, we put zeros only in the matrix. */
  Values[0] = 0.0;
  Values[1] = 0.0;
  Values[2] = 0.0;
  /* global ID's of local ID's */
  MyGlobalElements = Epetra_BlockMap_MyGlobalElements(bMap);

  /* At this point we simply set the nonzero structure of A.
   * Actual values will be inserted later (now all zeros)
   */
  for (i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0; 
      Indices[1] = 1;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      NumEntries = 3;
    }

    check = Epetra_CrsMatrix_InsertGlobalValues(A, MyGlobalElements[i],
        NumEntries, Values, Indices);
    assert(check == 0);
  }

  /* Finish up. */
  Epetra_CrsMatrix_FillComplete(A, TRUE);

  /*
   * ===========================================
   * E N D   O F   M A T R I X   C R E A T I O N
   * ===========================================
   */

  /*
   * Now the matrix STRUCTURE is set. We cannot add
   * new nonzero elements, but we can still change the
   * numerical values of all inserted elements (as we will
   * do later).
   */

  /*
   * =====================================================
   * B E G I N N I N G   O F  T H E   AM E S O S   P A R T
   * =====================================================
   */

  /*
   * For comments on the commands in this section, please
   * see file example_AmesosFactory.cpp.
   */

  Problem = Epetra_LinearProblem_Create();  

  rA = Epetra_RowMatrix_Cast(Epetra_CrsMatrix_Abstract(A));
  Epetra_LinearProblem_SetOperator_Matrix(Problem, rA);

  /* Initializes Amesos solver. Here we solve with Amesos_Klu. */

  Factory = Amesos_Create();
  Solver = Amesos_CreateSolver(Factory, "Amesos_Klu", Problem);

  /* Factory.Create() returns 0 if the requested solver
   * is not available
   */

/* 
  if (Solver == 0) {
    std::cerr << "Selected solver is not available" << std::endl;
    // return ok not to break the test harness
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }
*/

  /* At this point we can perform the numeric factorization.
   * Note that the matrix contains 0's only.
   */

  Amesos_BaseSolver_SymbolicFactorization(Solver);

  /* Now, we repopulate the matrix with entries corresponding
   * to a 1D Laplacian. LHS and RHS are still untouched.
   */

  for (i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0;   
      Indices[1] = 1;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements - 1;
      Indices[1] = NumGlobalElements - 2;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i] - 1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i] + 1;
      Values[0] = -1.0; 
      Values[1] = 2.0; 
      Values[2] = -1.0;
      NumEntries = 3;
    }

    check = Epetra_CrsMatrix_ReplaceGlobalValues(A, MyGlobalElements[i],
        NumEntries, Values, Indices);
    assert(check == 0);
  }

  /* ... and we can compute the numeric factorization. */
  Amesos_BaseSolver_NumericFactorization(Solver);

  /* Finally, we set up the LHS and the RHS vector (Random(). */
  b = Epetra_Vector_Create(bMap, TRUE);
  mb = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(b));
  Epetra_MultiVector_Random(mb);

  x = Epetra_Vector_Create(bMap, TRUE);
  mx = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(x));
  Epetra_MultiVector_PutScalar(mx, 0.0);

  Epetra_LinearProblem_SetLHS(Problem, mx);
  Epetra_LinearProblem_SetRHS(Problem, mb);
  
  Amesos_BaseSolver_Solve(Solver);

  /* Print out the timing information and get it from the solver */
  Amesos_BaseSolver_PrintTiming(Solver);
  
  TimingsList = Teuchos_ParameterList_Create();
  Amesos_BaseSolver_GetTiming(Solver, TimingsList);
  
  /* You can find out how much time was spent in ...
   * sfact_time, nfact_time, solve_time,
   * mtx_conv_time, mtx_redist_time, vec_redist_time */

  /* 1) The symbolic factorization
   *    (parameter doesn't always exist)
   */
  sfact_time = Teuchos_ParameterList_get_double_def( TimingsList, "Total symbolic factorization time", 0.0 );

  /* 2) The numeric factorization
   *    (always exists if NumericFactorization() is called)
   */
  nfact_time = Teuchos_ParameterList_get_double( TimingsList, "Total numeric factorization time" );

  /* 3) Solving the linear system
   *    (always exists if Solve() is called)
   */
  solve_time = Teuchos_ParameterList_get_double( TimingsList, "Total solve time" );

  /* 4) Converting the matrix to the accepted format for the solver
   *    (always exists if SymbolicFactorization() is called)
   */
  mtx_conv_time = Teuchos_ParameterList_get_double( TimingsList, "Total solve time" );

  /* 5) Redistributing the matrix for each solve to the accepted format for the solver */
  mtx_redist_time = Teuchos_ParameterList_get_double_def( TimingsList, "Total matrix redistribution time", 0.0 );

  /* 6) Redistributing the vector for each solve to the accepted format for the solver */
  vec_redist_time = Teuchos_ParameterList_get_double_def( TimingsList, "Total vector redistribution time", 0.0 );

  /*
   * ===========================================
   * E N D   O F   T H E   A M E S O S   P A R T
   * ===========================================
   */

  /*
   * ==================
   * compute ||Ax - b||
   * ==================
   */

  Ax = Epetra_Vector_Create(bMap, TRUE);
  mAx = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(Ax));
  Epetra_CrsMatrix_Multiply_Vector(A, FALSE, x, Ax);
  Epetra_MultiVector_Update_WithA(mAx, 1.0, mb, -1.0);
  Epetra_MultiVector_Norm2(mAx, &residual);

  if (!Epetra_Comm_MyPID(Comm))
    printf("After AMESOS solution, ||b-Ax||_2 = %g\n", residual);

  /* delete Solver. Do this before MPI_Finalize()
   * as MPI calls can occur in the destructor.
   */
  Amesos_BaseSolver_Destroy(&Solver);
    
  if (residual > 1e-5)
    return(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

#endif /* HAVE_CTRILINOS_AMESOS */

  return(EXIT_SUCCESS);

} /* end of main() */

