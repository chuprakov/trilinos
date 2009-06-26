#include "CEpetra_BlockMap.h"
#include "CEpetra_Distributor.h"
#include "Epetra_Export.h"
#include "CEpetra_Export.h"
#include "CEpetra_Export_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Cast )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  /* These casts should be allowed */
  ECHO(CT_Epetra_Object_ID_t objID = Epetra_Object_Cast(selfID));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, objID), true);
  ECHO(CT_Epetra_Export_ID_t impID = Epetra_Export_Cast(objID));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(objID, impID), true);
}

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Create )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Export_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Duplicate )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  /* Duplicate it */
  ECHO(CT_Epetra_Export_ID_t dupID = Epetra_Export_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_Export_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Destroy )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  ECHO(Epetra_Export_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , SourceMap )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID2 = Epetra_Export_SourceMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bsrcID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , TargetMap )
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

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(bsrcID, btarID));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t btarID2 = Epetra_Export_TargetMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(btarID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID );
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

