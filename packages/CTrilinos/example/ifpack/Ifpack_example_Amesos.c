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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "CTrilinos_enums.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_LinearProblem.h"
#include "CGaleri_Maps.h"
#include "CGaleri_CrsMatrices.h"
#include "CTeuchos_ParameterList.h"
#include "CAztecOO.h"
#include "CIfpack.h"

#include "az_aztec_defs.h"

#define CTRILINOS_CHK_ERR(ifpack_err) \
{ if (ifpack_err < 0) { \
  fprintf(stderr, "IFPACK ERROR %d, %s, line %d\n", \
    ifpack_err, __FILE__, __LINE__); \
    return(ifpack_err);  } }

int main(int argc, char *argv[])
{
  int NumProc, MyPID, nx, OverlapLevel;
  char PrecType[30];

  CT_Epetra_Comm_ID_t Comm;
  CT_Teuchos_ParameterList_ID_t GaleriList, List;
  CT_Epetra_Map_ID_t Map;
  CT_Epetra_Map_ID_t Map2;
  CT_Epetra_BlockMap_ID_t bMap2;
  CT_Epetra_CrsMatrix_ID_t A;
  CT_Epetra_RowMatrix_ID_t rA;
  CT_Ifpack_ID_t Factory;
  CT_Ifpack_Preconditioner_ID_t Prec;
  CT_Epetra_Vector_ID_t LHS, RHS;
  CT_Epetra_MultiVector_ID_t mLHS, mRHS;
  CT_Epetra_LinearProblem_ID_t Problem;
  CT_AztecOO_ID_t Solver;
  CT_Epetra_Operator_ID_t oPrec;

  /* initialize MPI and Epetra communicator */
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Comm = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(Epetra_MpiComm_Create(MPI_COMM_WORLD)));
#else
  Comm = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(Epetra_SerialComm_Create()));
#endif
  NumProc = Epetra_Comm_NumProc(Comm);
  MyPID = Epetra_Comm_MyPID(Comm);

  GaleriList = Teuchos_ParameterList_Create();
  
  /* The problem is defined on a 2D grid, global size is nx * nx. */
  nx = 30; 
  Teuchos_ParameterList_set_int(GaleriList, "nx", nx, "");
  Teuchos_ParameterList_set_int(GaleriList, "ny", nx * NumProc, "");
  Teuchos_ParameterList_set_int(GaleriList, "mx", 1, "");
  Teuchos_ParameterList_set_int(GaleriList, "my", NumProc, "");

  Map = Galeri_Maps_CreateMap("Cartesian2D", Comm, GaleriList);

  A = Galeri_CrsMatrices_CreateCrsMatrix("Laplace2D", Map, GaleriList);
  rA = Epetra_RowMatrix_Cast(Epetra_CrsMatrix_Abstract(A));

  /* ===============================================================
   * B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N
   * =============================================================== */

  List = Teuchos_ParameterList_Create();

  /* allocates an IFPACK factory. No data is associated 
   * to this object (only method Create()). */
  Factory = Ifpack_Create();

  /* create the preconditioner. For valid PrecType values,
   * please check the documentation */
  strcpy(PrecType, "Amesos");
  OverlapLevel = 2; /* must be >= 0. If Comm.NumProc() == 1, it is ignored. */

  Prec = Ifpack_CreatePreconditioner_UsingName(Factory, PrecType, rA, OverlapLevel);

  /* specify the Amesos solver to be used. 
   * If the selected solver is not available,
   * IFPACK will try to use Amesos' KLU (which is usually always
   * compiled). Amesos' serial solvers are:
   * "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu" */
  Teuchos_ParameterList_set_str(List, "amesos: solver type", "Amesos_Klu", "");

  /* sets the parameters */
  CTRILINOS_CHK_ERR(Ifpack_Preconditioner_SetParameters(Prec, List));

  /* initialize the preconditioner. At this point the matrix must
   * have been FillComplete()'d, but actual values are ignored.
   * At this call, Amesos will perform the symbolic factorization. */
  CTRILINOS_CHK_ERR(Ifpack_Preconditioner_Initialize(Prec));

  /* Builds the preconditioners, by looking for the values of 
   * the matrix. At this call, Amesos will perform the
   * numeric factorization. */
  CTRILINOS_CHK_ERR(Ifpack_Preconditioner_Compute(Prec));

  /* ===================================================
   * E N D   O F   I F P A C K   C O N S T R U C T I O N
   * =================================================== */

  /* At this point, we need some additional objects
   * to define and solve the linear system. */

  /* defines LHS and RHS */
  Map2 = Epetra_CrsMatrix_OperatorDomainMap(A);
  bMap2 = Epetra_BlockMap_Cast(Epetra_Map_Abstract(Map2));
  LHS = Epetra_Vector_Create(bMap2, FALSE);
  mLHS = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(LHS));
  RHS = Epetra_Vector_Create(bMap2, TRUE);
  mRHS = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(RHS));

  /* solution is constant */
  Epetra_MultiVector_PutScalar(mLHS, 1.0);
  /* now build corresponding RHS */
  Epetra_CrsMatrix_Apply(A, mLHS, mRHS);

  /* now randomize the solution */
  Epetra_MultiVector_Random(mRHS);

  /* need an Epetra_LinearProblem to define AztecOO solver */
  Problem = Epetra_LinearProblem_Create_FromMatrix(rA, mLHS, mRHS);

  /* now we can allocate the AztecOO solver */
  Solver = AztecOO_Create(Problem);

  /* specify solver */
  AztecOO_SetAztecOption(Solver, AZ_solver, AZ_gmres);
  AztecOO_SetAztecOption(Solver, AZ_output, 32);

  /* HERE WE SET THE IFPACK PRECONDITIONER */
  oPrec = Epetra_Operator_Cast(Ifpack_Preconditioner_Abstract(Prec));
  AztecOO_SetPrecOperator(Solver, oPrec);

  /* .. and here we solve
   * NOTE: with one process, the solver must converge in
   * one iteration. */
  AztecOO_Iterate_Current(Solver, 1550, 1e-8);

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
