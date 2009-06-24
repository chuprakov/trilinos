#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "Epetra_Vector.h"
#include "Epetra_DataAccess.h"
#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


/*
CT_Epetra_Comm_ID_t UnitTest_Create_Comm();
*/

TEUCHOS_UNIT_TEST( CEpetra_Helpers, CreateComm )
{
  ECHO(CT_Epetra_Comm_ID_t commID = UnitTest_Create_Comm());
}

