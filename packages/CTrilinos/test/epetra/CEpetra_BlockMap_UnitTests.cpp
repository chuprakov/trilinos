#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "Epetra_BlockMap.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 8);
  ECHO(int ElementSize = 2);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  /* This cast should be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID = Epetra_BlockMap_Cast(Epetra_BlockMap_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, bmapID), true);

  /* This cast should not be allowed */
  TEST_THROW(Epetra_Map_Cast(Epetra_BlockMap_Abstract(selfID)), Teuchos::m_bad_cast);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 9);
  ECHO(int ElementSize = 3);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_BlockMap_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);

  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(selfID), NumGlobalElements);
  TEST_EQUALITY(Epetra_BlockMap_ElementSize_Const(selfID), ElementSize);
  TEST_EQUALITY(Epetra_BlockMap_IndexBase(selfID), IndexBase);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Create_Linear )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumMyElements = 5);
  ECHO(int ElementSize = 2);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Linear(
       NumGlobalElements, NumMyElements, ElementSize, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_BlockMap_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Create_Arbitrary )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_BlockMap_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Create_Variable )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_BlockMap_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 8);
  ECHO(int ElementSize = 2);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(CT_Epetra_BlockMap_ID_t dupID = Epetra_BlockMap_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_BlockMap_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 8);
  ECHO(int ElementSize = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(Epetra_BlockMap_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MinAllGID )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {1+off, 2+off, 3+off, 0+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int val = Epetra_BlockMap_MinAllGID(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int val2 = 0);
  TEST_EQUALITY(val, val2);
}

/**********************************************************************
int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MaxAllGID )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {1+off, 2+off, 3+off, 0+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int val = Epetra_BlockMap_MaxAllGID(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int val2 = 3+NumMyElements*(NumProc-1));
  TEST_EQUALITY(val, val2);
}

/**********************************************************************
int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MinMyGID )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {1+off, 2+off, 3+off, 0+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int val = Epetra_BlockMap_MinMyGID(selfID));

  /* Now check the result of the call to the wrapper function */
  int val2 = MyGlobalElements[0];
  for (int i=1; i<NumMyElements; i++)
    if (val2 > MyGlobalElements[i]) val2 = MyGlobalElements[i];

  TEST_EQUALITY(val, val2);
}

/**********************************************************************
int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MaxMyGID )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {1+off, 2+off, 3+off, 0+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int val = Epetra_BlockMap_MaxMyGID(selfID));

  /* Now check the result of the call to the wrapper function */
  int val2 = MyGlobalElements[0];
  for (int i=1; i<NumMyElements; i++)
    if (val2 < MyGlobalElements[i]) val2 = MyGlobalElements[i];

  TEST_EQUALITY(val, val2);
}

/**********************************************************************
int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , NumGlobalElements )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 13);
  ECHO(int ElementSize = 3);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int els = Epetra_BlockMap_NumGlobalElements(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(els, NumGlobalElements);
}

/**********************************************************************
int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , NumMyElements )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  ECHO(int NumMyElements = MyPID);
  ECHO(int ElementSize = 1);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Linear(
       NumGlobalElements, NumMyElements, ElementSize, IndexBase, CommID));

  ECHO(int els = Epetra_BlockMap_NumMyElements(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(els, NumMyElements);
}

/**********************************************************************
int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , ElementSize_Const )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 12);
  ECHO(int ElementSize = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int sz = Epetra_BlockMap_ElementSize_Const(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(sz, ElementSize);
}

/**********************************************************************
int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , ElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  ECHO(int LID = 2);
  ECHO(int sz = Epetra_BlockMap_ElementSize(selfID, LID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(ElementSizeList[LID], sz);
}

/**********************************************************************
int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , IndexBase )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 9);
  ECHO(int ElementSize = 1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int ib = Epetra_BlockMap_IndexBase(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(ib, IndexBase);
}

/**********************************************************************
int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , NumGlobalPoints )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 13);
  ECHO(int ElementSize = 3);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int pts = Epetra_BlockMap_NumGlobalPoints(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int ngp = NumGlobalElements*ElementSize);
  TEST_EQUALITY(pts, ngp);
}

/**********************************************************************
int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , NumMyPoints )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 13);
  ECHO(int ElementSize = 3);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(int pts = Epetra_BlockMap_NumMyPoints(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int myels = Epetra_BlockMap_NumMyElements(selfID));

  ECHO(int nmp = myels*ElementSize);
  TEST_EQUALITY(pts, nmp);
}

/**********************************************************************
int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MinMyElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  ECHO(int sz = Epetra_BlockMap_MinMyElementSize(selfID));

  /* Now check the result of the call to the wrapper function */
  int min = ElementSizeList[0];
  for (int i=1; i<NumMyElements; i++)
    if (ElementSizeList[i] < min) min = ElementSizeList[i];

  TEST_EQUALITY(min, sz);
}

/**********************************************************************
int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MaxMyElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  ECHO(int sz = Epetra_BlockMap_MaxMyElementSize(selfID));

  /* Now check the result of the call to the wrapper function */
  int max = ElementSizeList[0];
  for (int i=1; i<NumMyElements; i++) /* sizes match on all processors */
    if (ElementSizeList[i] > max) max = ElementSizeList[i];

  TEST_EQUALITY(max, sz);
}

/**********************************************************************
int Epetra_BlockMap_MinElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MinElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  ECHO(int sz = Epetra_BlockMap_MinElementSize(selfID));

  /* Now check the result of the call to the wrapper function */
  int min = ElementSizeList[0];
  for (int i=1; i<NumMyElements; i++) /* sizes match on all processors */
    if (ElementSizeList[i] < min) min = ElementSizeList[i];

  TEST_EQUALITY(min, sz);
}

/**********************************************************************
int Epetra_BlockMap_MaxElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , MaxElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  int ElementSizeList[NumMyElements] = {2, 3, 2, 1};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Variable(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, CommID));

  ECHO(int sz = Epetra_BlockMap_MaxElementSize(selfID));

  /* Now check the result of the call to the wrapper function */
  int max = ElementSizeList[0];
  for (int i=1; i<NumMyElements; i++)
    if (ElementSizeList[i] > max) max = ElementSizeList[i];

  TEST_EQUALITY(max, sz);
}

/**********************************************************************
boolean Epetra_BlockMap_UniqueGIDs ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , UniqueGIDs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  int MyGlobalElements[NumMyElements] = {0, 1, 2, 3};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(boolean u = Epetra_BlockMap_UniqueGIDs(selfID));

  /* Now check the result of the call to the wrapper function */
  boolean u2 = FALSE;
  if (NumProc == 1) u2 = TRUE;
  TEST_EQUALITY(u, u2);
}

/**********************************************************************
boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , ConstantElementSize )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 12);
  ECHO(int ElementSize = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(boolean con = Epetra_BlockMap_ConstantElementSize(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(con, TRUE);
}

/**********************************************************************
boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , SameAs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  ECHO(int IndexBase = 0);

  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  int MyGlobalElements2[NumMyElements] = {0+off, 2+off, 1+off, 3+off};
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements2, ElementSize, IndexBase, CommID));

  ECHO(boolean same = Epetra_BlockMap_SameAs(selfID, MapID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(same, FALSE);
}

/**********************************************************************
boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , PointSameAs )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  ECHO(int IndexBase = 0);

  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  int MyGlobalElements2[NumMyElements] = {0+off, 2+off, 1+off, 3+off};
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements2, ElementSize, IndexBase, CommID));

  ECHO(boolean same = Epetra_BlockMap_PointSameAs(selfID, MapID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(same, TRUE);
}

/**********************************************************************
boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , LinearMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(boolean lin = Epetra_BlockMap_LinearMap(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(lin, FALSE);
}

/**********************************************************************
boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , DistributedGlobal )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(boolean u = Epetra_BlockMap_DistributedGlobal(selfID));

  /* Now check the result of the call to the wrapper function */
  boolean u2 = TRUE;
  if (NumProc == 1) u2 = FALSE;
  TEST_EQUALITY(u, u2);
}

/**********************************************************************
int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_IsOneToOne ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , IsOneToOne )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  int MyGlobalElements[NumMyElements] = {0, 1, 2, 3};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID));

  ECHO(boolean u = Epetra_BlockMap_IsOneToOne(selfID));

  /* Now check the result of the call to the wrapper function */
  boolean u2 = FALSE;
  if (NumProc == 1) u2 = TRUE;
  TEST_EQUALITY(u, u2);
}

/**********************************************************************
void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_BlockMap , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create a map to duplicate */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int IndexBase = 0);
  ECHO(int ElementSize1 = 3);
  ECHO(int NumGlobalElements1 = 4);
  ECHO(CT_Epetra_BlockMap_ID_t  mapID = Epetra_BlockMap_Create(
       NumGlobalElements1, ElementSize1, IndexBase, CommID));

  /* Create the one to operate on */
  ECHO(int ElementSize2 = 2);
  ECHO(int NumGlobalElements2 = 6);
  ECHO(CT_Epetra_BlockMap_ID_t selfID = Epetra_BlockMap_Create(
       NumGlobalElements2, ElementSize2, IndexBase, CommID));

  /* Check the initial state */
  ECHO(int els = Epetra_BlockMap_NumGlobalElements(selfID));
  TEST_EQUALITY(els, NumGlobalElements2);

  /* Test out the wrapper and check that it worked */
  ECHO(Epetra_BlockMap_Assign(selfID, mapID));
  ECHO(els = Epetra_BlockMap_NumGlobalElements(selfID));
  TEST_EQUALITY(els, NumGlobalElements1);
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

