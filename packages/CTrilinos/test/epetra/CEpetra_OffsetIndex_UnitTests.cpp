#include "CEpetra_CrsGraph.h"
#include "CEpetra_Export.h"
#include "CEpetra_Import.h"
#include "Epetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
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

