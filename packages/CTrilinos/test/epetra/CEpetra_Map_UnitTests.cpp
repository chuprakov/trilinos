#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 13);
  ECHO(int ElementSize = 4);
  ECHO(int IndexBase = 0);

  /* This cast should be allowed */
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_Map_ID_t dupID = Epetra_Map_Cast(
       selfID));

  /* This cast should be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID = Epetra_BlockMap_Cast(
       selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bmapID), NumGlobalElements);
  ECHO(CT_Epetra_Map_ID_t mapID = Epetra_Map_Cast(
       bmapID));

  /* This cast should not be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID2 = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));
  TEST_THROW(Epetra_Map_Cast(bmapID2), Teuchos::m_bad_cast);
}


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create_Linear )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumMyElements = 5);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create_Arbitrary )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[] = {0+off, 1+off, 2+off, 3+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int NumGlobalElements = 4);
  CT_Epetra_Map_ID_t selfID = UnitTest_Create_Map(NumGlobalElements);

  ECHO(CT_Epetra_Map_ID_t dupID = Epetra_Map_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
}


/**********************************************************************
void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int NumGlobalElements = 6);
  CT_Epetra_Map_ID_t selfID = UnitTest_Create_Map(NumGlobalElements);

  ECHO(Epetra_Map_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
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


