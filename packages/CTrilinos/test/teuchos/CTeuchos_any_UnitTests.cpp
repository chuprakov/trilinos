#include "CTrilinos_config.h"
#include "Teuchos_any.hpp"
#include "CTeuchos_any.h"
#include "CTeuchos_any_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {



/**********************************************************************
CT_Teuchos_any_ID_t Teuchos_any_Create (  );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID );
 **********************************************************************/

/**********************************************************************
void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );
 **********************************************************************/

/**********************************************************************
void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );
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

