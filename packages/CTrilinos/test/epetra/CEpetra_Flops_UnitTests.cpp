#include "Epetra_Flops.h"
#include "CEpetra_Flops.h"
#include "CEpetra_Flops_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_Flops_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );
 **********************************************************************/

/**********************************************************************
double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );
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

