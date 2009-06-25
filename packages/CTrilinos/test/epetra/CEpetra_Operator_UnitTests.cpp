#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "CEpetra_Operator.h"
#include "CEpetra_Operator_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Operator_ID_t Epetra_Operator_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );
 **********************************************************************/

/**********************************************************************
int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
const char * Epetra_Operator_Label ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_Operator_HasNormInf ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );
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

