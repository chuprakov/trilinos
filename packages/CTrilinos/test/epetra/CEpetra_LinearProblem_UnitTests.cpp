#include "CTrilinos_config.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Operator.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "CEpetra_LinearProblem.h"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {



/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CTrilinos_Object_ID_t Epetra_LinearProblem_Abstract ( 
  CT_Epetra_LinearProblem_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );
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

