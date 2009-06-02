#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "Epetra_DataAccess.h"
#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_BlockMap_ID_t MapID = UnitTest_Create_BlockMapCast(9);

  /* This cast should be allowed */
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create_AndZero(
       MapID));
  ECHO(CT_Epetra_Vector_ID_t dupID = Epetra_Vector_Cast(selfID));

  /* This cast should be allowed */
  ECHO(CT_Epetra_MultiVector_ID_t mvecID = Epetra_MultiVector_Cast(selfID));
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Cast(mvecID));

  /* This cast should not be allowed */
  ECHO(int NumVectors = 3);
  ECHO(CT_Epetra_MultiVector_ID_t mvecID2 = Epetra_MultiVector_Create_AndZero(
       MapID, NumVectors));
  TEST_THROW(Epetra_Vector_Cast(mvecID2), Teuchos::m_bad_cast);

  /* Now check the result of the call to the wrapper function */
  /* ... */
}


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_BlockMap_ID_t MapID = UnitTest_Create_BlockMapCast(7);

  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(
       MapID, zeroOut));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Create_AndZero ( 
  CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create_AndZero )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_BlockMap_ID_t MapID = UnitTest_Create_BlockMapCast(9);

  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create_AndZero(
       MapID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_Vector_ID_t SourceID = UnitTest_Create_Vector();

  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Duplicate(
       SourceID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 1);
}


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * V );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create_FromArray )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_BlockMap_ID_t MapID = UnitTest_Create_BlockMapCast(6);
  double V[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create_FromArray(
       Copy, MapID, V));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , FromSource )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_BlockMap_ID_t MapID = UnitTest_Create_BlockMapCast(9);

  ECHO(int NumVectors = 3);
  ECHO(CT_Epetra_MultiVector_ID_t SourceID = Epetra_MultiVector_Create_AndZero(
       MapID, NumVectors));

  ECHO(int Index = 1);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_FromSource(
       View, SourceID, Index));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}


/**********************************************************************
void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  CT_Epetra_Vector_ID_t selfID = UnitTest_Create_Vector();

  ECHO(Epetra_Vector_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}


/**********************************************************************
int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V );
 **********************************************************************/


/**********************************************************************
int Epetra_Vector_ResetView ( 
  CT_Epetra_Vector_ID_t selfID, double * Values_in );
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


