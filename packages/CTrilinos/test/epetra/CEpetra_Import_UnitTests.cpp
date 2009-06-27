#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CEpetra_Distributor.h"
#include "Epetra_Import.h"
#include "CEpetra_Import.h"
#include "CEpetra_Import_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  /* These casts should be allowed */
  ECHO(CT_Epetra_Object_ID_t objID = Epetra_Object_Cast(selfID));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, objID), true);
  ECHO(CT_Epetra_Import_ID_t impID = Epetra_Import_Cast(objID));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(objID, impID), true);
}

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Import_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
  ECHO(Teuchos::RCP<Epetra_Import> rcp1 = CEpetra::getImport(selfID));
  TEST_EQUALITY_CONST(rcp1.is_null(), false);
}

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  /* Duplicate it */
  ECHO(CT_Epetra_Import_ID_t dupID = Epetra_Import_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_Import_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  ECHO(Epetra_Import_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumSameIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumSameIDs = Epetra_Import_NumSameIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumSameIDs, 3);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumSameIDs, 0);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumSameIDs, 0);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumPermuteIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumPermuteIDs = Epetra_Import_NumPermuteIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumPermuteIDs, 0);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , PermuteFromLIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumPermuteIDs = Epetra_Import_NumPermuteIDs(selfID));
    ECHO(int *ret = Epetra_Import_PermuteFromLIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumPermuteIDs, 0);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      TEST_EQUALITY_CONST(ret[0], 1);
      TEST_EQUALITY_CONST(ret[1], 2);
      TEST_EQUALITY_CONST(ret[2], 3);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      TEST_EQUALITY_CONST(ret[0], 2);
      TEST_EQUALITY_CONST(ret[1], 3);
      TEST_EQUALITY_CONST(ret[2], 4);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , PermuteToLIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumPermuteIDs = Epetra_Import_NumPermuteIDs(selfID));
    ECHO(int *ret = Epetra_Import_PermuteToLIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumPermuteIDs, 0);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 1);
      TEST_EQUALITY_CONST(ret[2], 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumPermuteIDs, 3);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 1);
      TEST_EQUALITY_CONST(ret[2], 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumRemoteIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumRemoteIDs = Epetra_Import_NumRemoteIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , RemoteLIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumRemoteIDs = Epetra_Import_NumRemoteIDs(selfID));
    ECHO(int *ret = Epetra_Import_RemoteLIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 3);
      TEST_EQUALITY_CONST(ret[1], 4);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 4);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumRemoteIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 1);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumExportIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumExportIDs = Epetra_Import_NumExportIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , ExportLIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumExportIDs = Epetra_Import_NumExportIDs(selfID));
    ECHO(int *ret = Epetra_Import_ExportLIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , ExportPIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumExportIDs = Epetra_Import_NumExportIDs(selfID));
    ECHO(int *ret = Epetra_Import_ExportPIDs(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 1);
      TEST_EQUALITY_CONST(ret[1], 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumExportIDs, 2);
      TEST_EQUALITY_CONST(ret[0], 0);
      TEST_EQUALITY_CONST(ret[1], 1);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumSend )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumSend = Epetra_Import_NumSend(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumSend, 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumSend, 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumSend, 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

#ifdef HAVE_MPI

TEUCHOS_UNIT_TEST ( Epetra_Import , NumRecv )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  TEST_EQUALITY_CONST(NumProc, 3);

  if (NumProc == 3) {
    /* Create the source map */
    ECHO(int IndexBase = 0);
    ECHO(const int NumMyElements = 3);
    ECHO(int NumGlobalElements = NumMyElements * NumProc);
    ECHO(int off = NumMyElements*MyPID);
    int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
    ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

    /* Create the target map */
    ECHO(const int NumMyElements2 = 5);
    ECHO(int NumGlobalElements2 = NumMyElements2 * NumProc);
    int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
    int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
    int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
    int MyGlobalElements2[NumMyElements2];
    for (int i=0; i<NumProc; i++) {
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
    ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
         NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID));
    ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

    /* Create an importer */
    ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

    /* Try out the wrapper function */
    ECHO(int NumRecv = Epetra_Import_NumRecv(selfID));

    /* Now check the result of the call to the wrapper function */
    switch (MyPID) {
    case 0:
      TEST_EQUALITY_CONST(NumRecv, 2);
      break;
    case 1:
      TEST_EQUALITY_CONST(NumRecv, 2);
      break;
    case 2:
      TEST_EQUALITY_CONST(NumRecv, 2);
      break;
    }
  }
}

#endif /* HAVE_MPI */

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , SourceMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID2 = Epetra_Import_SourceMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bsrcID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , TargetMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t btarID2 = Epetra_Import_TargetMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(btarID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Distributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Cast(srcID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Cast(tarID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID));

  ECHO(CT_Epetra_Distributor_ID_t dID = Epetra_Import_Distributor(selfID));


  /* ??? WHY DOESN'T THIS WORK ??? */


  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dID.type, CT_Epetra_Distributor_ID);
  TEST_EQUALITY_CONST(dID.index, 0);
  TEST_INEQUALITY(CEpetra::getDistributor(dID), Teuchos::null);
}


/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

