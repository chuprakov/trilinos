#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )


#ifdef EPETRA_MPI
#define COMM_NAME Epetra_MpiComm
#include "CEpetra_MpiComm.h"
#else /* EPETRA_MPI */
#define COMM_NAME Epetra_SerialComm
#include "CEpetra_SerialComm.h"
#endif /* EPETRA_MPI */
#define COMM_TYPE CLASS_TYPE(COMM_NAME)
#define COMM_ENUM CLASS_ENUM(COMM_NAME)
#define COMM_CREATE BUILD_CALL(COMM_NAME, Create)


namespace {


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

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
  ECHO(CT_Epetra_Map_ID_t mapID = Epetra_Map_Cast(
       bmapID));

  /* This cast should not be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID2 = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));
  TEST_THROW(CT_Epetra_Map_ID_t mapID2 = Epetra_Map_Cast(
       bmapID2), Teuchos::m_bad_cast);

  /* Now check the result of the call to the wrapper function */
  /* ... */
}


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

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
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

  ECHO(int NumMyElements = 5);
  ECHO(int NumProc = BUILD_CALL(COMM_NAME, NumProc) (cid));
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
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

  ECHO(int NumMyElements = 4);
  ECHO(int NumProc = BUILD_CALL(COMM_NAME, NumProc) (cid));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = BUILD_CALL(COMM_NAME, MyPID) (cid));
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

  /* Create everything we need to pass to the constructor */
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));

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

  /* Create everything we need to pass to the constructor */
  ECHO(COMM_TYPE cid = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(cid));

  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));

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


