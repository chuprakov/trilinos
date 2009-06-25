#include "CEpetra_Flops.h"
#include "CEpetra_Comm.h"
#include "Epetra_CompObject.h"
#include "CEpetra_CompObject.h"
#include "CEpetra_CompObject_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_CompObject_ID_t Epetra_CompObject_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CompObject , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t selfID = Epetra_CompObject_Create());

  /* This cast should be allowed */
  ECHO(CT_Epetra_CompObject_ID_t coID = Epetra_CompObject_Cast(selfID));

  /* This cast should not be allowed */
  TEST_THROW(Epetra_Comm_Cast(selfID), Teuchos::m_bad_cast);
}

/**********************************************************************
CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CompObject , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t selfID = Epetra_CompObject_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_CompObject_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CompObject , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t selfID = Epetra_CompObject_Create());
  ECHO(CT_Epetra_CompObject_ID_t dupID = Epetra_CompObject_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_CompObject_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
}

/**********************************************************************
void Epetra_CompObject_Destroy ( 
  CT_Epetra_CompObject_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CompObject , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t selfID = Epetra_CompObject_Create());
  ECHO(Epetra_CompObject_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in );
 **********************************************************************/

/**********************************************************************
void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID );
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

