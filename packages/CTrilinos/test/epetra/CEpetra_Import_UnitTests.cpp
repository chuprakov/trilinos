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

/**********************************************************************
int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

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

