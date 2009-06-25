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

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID );
 **********************************************************************/

/**********************************************************************
void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );
 **********************************************************************/

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

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

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

