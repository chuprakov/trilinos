#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AMESOS

#include "CAmesos_BaseSolver.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"
#include "Amesos.h"
#include "CAmesos.h"
#include "CAmesos_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {



/**********************************************************************
CT_Amesos_ID_t Amesos_Create (  );
 **********************************************************************/

/**********************************************************************
void Amesos_Destroy ( CT_Amesos_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( 
  CT_Amesos_ID_t selfID, const char * ClassType, 
  CT_Epetra_LinearProblem_ID_t LinearProblemID );
 **********************************************************************/

/**********************************************************************
boolean Amesos_Query ( 
  CT_Amesos_ID_t selfID, const char * ClassType );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );
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

#endif /* HAVE_CTRILINOS_AMESOS */

