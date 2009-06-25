#include "Epetra_Object.h"
#include "CEpetra_Object.h"
#include "CEpetra_Object_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );
 **********************************************************************/

/**********************************************************************
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );
 **********************************************************************/

/**********************************************************************
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );
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

