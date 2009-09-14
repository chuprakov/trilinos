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
#include "CEpetra_Comm.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#endif
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CEpetra_Import.h"
#include "CEpetra_Export.h"
#include "CTrilinos_enums.h"


int checkresult(int MyResult, int Correct, CT_Epetra_Comm_ID_t CommID, int MyPID, const char *pc)
{
  int Failed = (MyResult == Correct ? 0 : 1);

  if (Failed != 0) {
    printf("  %s FAILED on processor %d:  %d != %d\n", pc, MyPID, MyResult, Correct);
    fflush(stdout);
  }

  int SomeFailed;
  Epetra_Comm_MaxAll_Int(CommID, &Failed, &SomeFailed, 1);

  if ((SomeFailed == 0) && (MyPID == 0)) {
    printf("  %s passed!\n", pc, MyPID);
    fflush(stdout);
  }
  Epetra_Comm_Barrier(CommID);

  return SomeFailed;
}

int checkresultarray(int MyResult[], int Correct[], int Count, CT_Epetra_Comm_ID_t CommID, int MyPID, const char *pc)
{
  int Failed = 0;
  for (int i=0; i<Count; i++) {
    if (MyResult[i] != Correct[i]) Failed = 1;
  }

  if (Failed != 0) {
    printf("  %s FAILED on processor %d: {%d", pc, MyPID, MyResult[0]);
    for (int i=1; i<Count; i++) printf(",%d", MyResult[i]);
    printf("} != {%d", Correct[0]);
    for (int i=1; i<Count; i++) printf(",%d", Correct[i]);
    printf("}\n");
    fflush(stdout);
  }

  int SomeFailed;
  Epetra_Comm_MaxAll_Int(CommID, &Failed, &SomeFailed, 1);

  if ((SomeFailed == 0) && (MyPID == 0)) {
    printf("  %s passed!\n", pc, MyPID);
    fflush(stdout);
  }
  Epetra_Comm_Barrier(CommID);

  return SomeFailed;
}

int main(int argc, char *argv[])
{
  int success = 0;

#ifdef HAVE_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);

  const int ForceNumProc = 3;

  /* Set up communication */
  CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(
      Epetra_MpiComm_Create( MPI_COMM_WORLD )));

  int MyPID = Epetra_Comm_MyPID(CommID);
  int NumProc = Epetra_Comm_NumProc(CommID);

  if (NumProc != ForceNumProc) {
    fprintf(stderr, "Run this with %d processors!\n", ForceNumProc);
    printf("End Result: TEST FAILED\n");
    success = 1;
    return success;
  }

  Epetra_Comm_Barrier(CommID);

  /* Create the source map */
  if (MyPID == 0) printf("Creating source map...\n");
  int IndexBase = 0;
  const int NumMyElements = 3;
  int NumGlobalElements = NumMyElements * NumProc;
  int off = NumMyElements*MyPID;
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
  CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID);
  CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(srcID));

  Epetra_Comm_Barrier(CommID);

  /* Create the target map */
  if (MyPID == 0) printf("Creating target map...\n");
  const int NumMyElements2 = 5;
  int NumGlobalElements2 = NumMyElements2 * NumProc;
  int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
  int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
  int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
  int MyGlobalElements2[NumMyElements2];
  for (int i=0; i<NumMyElements2; i++) {
    switch (MyPID) {
    case 0:
      MyGlobalElements2[i] = MyGlobalElements2a[i];
      break;
    case 1:
      MyGlobalElements2[i] = MyGlobalElements2b[i];
      break;
    case 2:
      MyGlobalElements2[i] = MyGlobalElements2c[i];
      break;
    }
  }
  CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID);
  CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(tarID));

  Epetra_Comm_Barrier(CommID);

{ /* IMPORT BLOCK */

  /* Create an importer */
  if (MyPID == 0) printf("Creating importer...\n");
  CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID);

  Epetra_Comm_Barrier(CommID);

  /* Try out the wrapper functions... */
  if (MyPID == 0) printf("Testing importer...\n");

  Epetra_Comm_Barrier(CommID);

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumSameIDs(selfID);
    int Correct[ForceNumProc] = {3, 0, 0};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSameIDs");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumPermuteIDs(selfID);
    int Correct[ForceNumProc] = {0, 3, 3};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumPermuteIDs");
    success += success1;

    /* The Trilinos 9.0 Doxygen page example for Epetra_Import has the results of
     * PermuteToLIDs and PermuteFromLIDs switched! */

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_PermuteToLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {1, 2, 3};
      int Correct2c[MaxCount] = {2, 3, 4};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteToLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_PermuteFromLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {0, 1, 2};
      int Correct2c[MaxCount] = {0, 1, 2};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteFromLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumRemoteIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRemoteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_RemoteLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {3, 4};
      int Correct2b[MaxCount] = {0, 4};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "RemoteLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumExportIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumExportIDs");
    success += success1;

    /* The Trilinos 9.0 Doxygen page example for Epetra_Import has the results of
     * ExportLIDs in the wrong order for processors 0 and 2! */

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_ExportLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {2, 0};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {2, 0};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_ExportPIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {1, 2};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportPIDs");
    }
  }

  /****************************************************************/

  {
    /* The source code for Epetra_Import never sets NumSend to non-zero! */

    int MyResult = Epetra_Import_NumSend(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSend");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumRecv(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRecv");
  }

  /****************************************************************/

  Epetra_Comm_Barrier(CommID);

  /* Destroy importer */
  Epetra_Import_Destroy(&selfID);

  Epetra_Comm_Barrier(CommID);

} /* IMPORT BLOCK */

{ /* EXPORT BLOCK */

  /* Create an exporter */
  if (MyPID == 0) printf("Creating exporter...\n");
  /* Intentionally swapping source and target maps! */
  CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(btarID, bsrcID);

  Epetra_Comm_Barrier(CommID);

  /* Try out the wrapper functions... */
  if (MyPID == 0) printf("Testing exporter...\n");

  Epetra_Comm_Barrier(CommID);

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumSameIDs(selfID);
    int Correct[ForceNumProc] = {3, 0, 0};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSameIDs");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumPermuteIDs(selfID);
    int Correct[ForceNumProc] = {0, 3, 3};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumPermuteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_PermuteToLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {0, 1, 2};
      int Correct2c[MaxCount] = {0, 1, 2};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteToLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_PermuteFromLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {1, 2, 3};
      int Correct2c[MaxCount] = {2, 3, 4};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteFromLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumRemoteIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRemoteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_RemoteLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {2, 0};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {2, 0};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "RemoteLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumExportIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumExportIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_ExportLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {3, 4};
      int Correct2b[MaxCount] = {0, 4};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_ExportPIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {1, 2};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportPIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumSend(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSend");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumRecv(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRecv");
  }

  /****************************************************************/

  Epetra_Comm_Barrier(CommID);

  /* Destroy exporter */
  Epetra_Export_Destroy(&selfID);

  Epetra_Comm_Barrier(CommID);

} /* EXPORT BLOCK */

  if (success == 0)
    printf( "\nEnd Result: TEST PASSED\n" );
  else
    printf( "\nEnd Result: TEST FAILED\n" );
  
  MPI_Finalize() ;

#else /* HAVE_MPI */

  printf("End Result: TEST FAILED\n");
  return 1;

#endif /* HAVE_MPI */

  return ((success == 0) ? 0 : 1);
}

