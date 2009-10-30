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

#include <stdio.h>

#include "az_aztec_defs.h"
#include "CAztecOO.h"

#include "CTrilinos_enums.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "CEpetra_MpiComm.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"

#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "CEpetra_CrsMatrix.h"

int main(int argc, char *argv[])
{
  CT_Epetra_Comm_ID_t Comm;
  CT_Epetra_Map_ID_t Map;
  CT_Epetra_BlockMap_ID_t bMap;
  CT_Epetra_CrsMatrix_ID_t A;
  CT_Epetra_Vector_ID_t x, b;
  CT_Epetra_MultiVector_ID_t mx, mb;
  CT_Epetra_RowMatrix_ID_t rA;
  CT_Epetra_LinearProblem_ID_t problem;
  CT_AztecOO_ID_t solver;

  int NumMyElements, NumGlobalElements, i, GlobalRow, RowLess1, RowPlus1;
  double negOne, posTwo;

#ifdef HAVE_MPI
  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  Comm = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(Epetra_MpiComm_Create(MPI_COMM_WORLD)));
#else
  Comm = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(Epetra_SerialComm_Create()));
#endif

  NumMyElements = 1000;

  /* Construct a Map that puts same number of equations on each processor */
  Map = Epetra_Map_Create_Linear(-1, NumMyElements, 0, Comm);
  bMap = Epetra_BlockMap_Cast(Epetra_Map_Abstract(Map));

  NumGlobalElements = Epetra_BlockMap_NumGlobalElements(bMap);

  /* Create a Epetra_Matrix */
  A = Epetra_CrsMatrix_Create(CT_Epetra_DataAccess_E_Copy, Map, 3, FALSE);
  
  /* Add rows one-at-a-time */
  negOne = -1.0;
  posTwo = 2.0;
  for (i=0; i<NumMyElements; i++) {
    GlobalRow = Epetra_CrsMatrix_GRID(A, i);
    RowLess1 = GlobalRow - 1;
    RowPlus1 = GlobalRow + 1;

    if (RowLess1 != -1)
      Epetra_CrsMatrix_InsertGlobalValues(A, GlobalRow, 1, &negOne, &RowLess1);
    if (RowPlus1 != NumGlobalElements)
      Epetra_CrsMatrix_InsertGlobalValues(A, GlobalRow, 1, &negOne, &RowPlus1);
    Epetra_CrsMatrix_InsertGlobalValues(A, GlobalRow, 1, &posTwo, &GlobalRow);
  }
  
  /* Finish up */
  Epetra_CrsMatrix_FillComplete(A, TRUE);

  /* Create x and b vectors */
  x = Epetra_Vector_Create(bMap, TRUE);
  b = Epetra_Vector_Create(bMap, TRUE);

  /* Cast vectors to multivectors */
  mx = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(x));
  mb = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(b));
  Epetra_MultiVector_Random(mb);

  /* Cast from CrsMatrix to RowMatrix */
  rA = Epetra_RowMatrix_Cast(Epetra_CrsMatrix_Abstract(A));

  /* Create Linear Problem */
  problem = Epetra_LinearProblem_Create_FromMatrix(rA, mx, mb);

  /* Create AztecOO instance */
  solver = AztecOO_Create_FromLinearProblem(problem);

  AztecOO_SetAztecOption(solver, AZ_precond, AZ_Jacobi);
  AztecOO_Iterate_Current(solver, 1000, 1.0E-8);

  printf("Solver performed %d iterations.\n", AztecOO_NumIters(solver));
  printf("Norm of true residual = %g\n", AztecOO_TrueResidual(solver));

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return 0;
}

