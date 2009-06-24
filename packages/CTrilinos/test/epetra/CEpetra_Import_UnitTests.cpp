#include "CEpetra_BlockMap.h"
#include "CEpetra_Distributor.h"
#include "Epetra_Import.h"
#include "CEpetra_Import.h"
#include "CEpetra_Import_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

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

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID );
 **********************************************************************/

/**********************************************************************
void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );
 **********************************************************************/

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

