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
#include "CTrilinos_utils_templ.hpp"

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

TEUCHOS_UNIT_TEST( Epetra_Operator , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create the source matrix */
  ECHO(int NumIndicesPerRow = 4);
  ECHO(Epetra_DataAccess CV = Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t crsID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Cast it to an operator */
  ECHO(CT_Epetra_Operator_ID_t selfID = Epetra_Operator_Cast(Epetra_CrsMatrix_Abstract(crsID)));

  ECHO(Epetra_Operator_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , SetUseTranspose )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create the source matrix */
  ECHO(int NumIndicesPerRow = 4);
  ECHO(Epetra_DataAccess CV = Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t crsID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Cast it to an operator */
  ECHO(CT_Epetra_Operator_ID_t selfID = Epetra_Operator_Cast(Epetra_CrsMatrix_Abstract(crsID)));

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(selfID, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(selfID));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(selfID, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(selfID));
  TEST_EQUALITY(tr, tr2);
}

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

TEUCHOS_UNIT_TEST( Epetra_Operator , UseTranspose )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create the source matrix */
  ECHO(int NumIndicesPerRow = 4);
  ECHO(Epetra_DataAccess CV = Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t crsID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Cast it to an operator */
  ECHO(CT_Epetra_Operator_ID_t selfID = Epetra_Operator_Cast(Epetra_CrsMatrix_Abstract(crsID)));

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(selfID, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(selfID));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(selfID, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(selfID));
  TEST_EQUALITY(tr, tr2);
}

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

