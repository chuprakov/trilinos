#include "CEpetra_Comm.h"
#include "Epetra_Time.h"
#include "CEpetra_Time.h"
#include "CEpetra_Time_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {



/**********************************************************************
CT_Epetra_Time_ID_t Epetra_Time_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Time_ID_t Epetra_Time_Create ( 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( 
  CT_Epetra_Time_ID_t TimeID );
 **********************************************************************/

/**********************************************************************
double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_Time_Assign ( 
  CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID );
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

