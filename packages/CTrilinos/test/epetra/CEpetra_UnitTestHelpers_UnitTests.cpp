#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "Epetra_DataAccess.h"
#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


/*
CT_Epetra_Comm_ID_t UnitTest_Create_Comm();
*/

TEUCHOS_UNIT_TEST( CEpetra_Helpers, CreateCommCast )
{
  ECHO(CT_Epetra_Comm_ID_t commID = UnitTest_Create_Comm());
}

/*
CT_Epetra_Map_ID_t UnitTest_Create_Map(int NumGlobalElements);
*/

TEUCHOS_UNIT_TEST( CEpetra_Helpers, CreateMap )
{
  ECHO(CT_Epetra_Map_ID_t mapID = UnitTest_Create_Map(7));
}

/*
CT_Epetra_BlockMap_ID_t UnitTest_Create_BlockMapCast(int NumGlobalElements);
*/

TEUCHOS_UNIT_TEST( CEpetra_Helpers, CreateBlockMapCast )
{
  ECHO(CT_Epetra_BlockMap_ID_t mapID = UnitTest_Create_BlockMapCast(8));
}

/*
CT_Epetra_Vector_ID_t UnitTest_Create_Vector(int NumGlobalElements = 7);
*/

TEUCHOS_UNIT_TEST( CEpetra_Helpers, CreateVector )
{
  ECHO(CT_Epetra_Vector_ID_t vecID = UnitTest_Create_Vector());
}

