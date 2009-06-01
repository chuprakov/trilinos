#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )


#ifdef EPETRA_MPI
#define COMM_NAME Epetra_MpiComm
#include "CEpetra_MpiComm.h"
#else /* EPETRA_MPI */
#define COMM_NAME Epetra_SerialComm
#include "CEpetra_SerialComm.h"
#endif /* EPETRA_MPI */
#define COMM_TYPE CLASS_TYPE(COMM_NAME)
#define COMM_ENUM CLASS_ENUM(COMM_NAME)
#define COMM_CREATE BUILD_CALL(COMM_NAME, Create)


namespace {


/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_Comm_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Comm , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(COMM_TYPE id = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t selfID = Epetra_Comm_Cast(id));

  /* If no exception thrown, then test was successful */
}


/**********************************************************************
void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Comm , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(COMM_TYPE id = COMM_CREATE ());
  ECHO(CT_Epetra_Comm_ID_t selfID = Epetra_Comm_Cast(id));

  ECHO(Epetra_Comm_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);

  /* Also try one that should fail since wrong type */
  TEST_THROW(Epetra_Comm_Destroy(&id), CTrilinos::CTrilinosTypeMismatchError);
}


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


