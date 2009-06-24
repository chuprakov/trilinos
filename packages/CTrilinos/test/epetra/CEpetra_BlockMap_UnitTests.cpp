#include "CEpetra_BlockMapData.h"
#include "CEpetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

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
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID );
 **********************************************************************/

/**********************************************************************
void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID );
 **********************************************************************/

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

/**********************************************************************
int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

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

/**********************************************************************
int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MinElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_BlockMap_MaxElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_UniqueGIDs ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID );
 **********************************************************************/

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

/**********************************************************************
void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMapData_ID_t Epetra_BlockMap_DataPtr ( 
  CT_Epetra_BlockMap_ID_t selfID );
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

