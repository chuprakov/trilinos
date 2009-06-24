#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_MultiVector.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Vector.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Vector_Cpp.hpp"
#include "Teuchos_RCP.hpp"
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
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(MapID, zeroOut));

  /* These casts should be allowed */
  ECHO(CT_Epetra_Vector_ID_t dupID = Epetra_Vector_Cast(selfID));
  ECHO(CT_Epetra_MultiVector_ID_t mvecID = Epetra_MultiVector_Cast(selfID));
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Cast(mvecID));

  /* This cast should not be allowed */
  ECHO(int NumVectors = 3);
  ECHO(CT_Epetra_MultiVector_ID_t mvecID2 = Epetra_MultiVector_Create(
       MapID, NumVectors, TRUE));
  TEST_THROW(Epetra_Vector_Cast(mvecID2), Teuchos::m_bad_cast);
}

/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(MapID, zeroOut));

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
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 8);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_Vector_ID_t SourceID = Epetra_Vector_Create(MapID, zeroOut));

  /* Call the wrapper function */
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Duplicate(SourceID));

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
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  double V[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create_FromArray(Copy, MapID, V));

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
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  ECHO(int NumVectors = 3);
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_MultiVector_ID_t SourceID = Epetra_MultiVector_Create(MapID, NumVectors, zeroOut));
  ECHO(int Index = 1);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_FromSource(View, SourceID, Index));

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
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 13);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(MapID, zeroOut));

  /* Call the wrapper function */
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

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceGlobalValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 11);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(int NumEntries = 5);
  double vals[] = {2.0, 7.4, 21.0, 8.5, 6.7};
  int inds[] = {0, 1, 4, 8, 10};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_ReplaceGlobalValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->ReplaceGlobalValues(NumEntries, vals, inds));

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceMyValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 5);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t mID = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(mID));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(int NumEntries = 3);
  double vals[] = {1.0, 5.2, 2.1};
  int inds[] = {0, 1, 4};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_ReplaceMyValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->ReplaceMyValues(NumEntries, vals, inds));

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoGlobalValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 11);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID)));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(int NumEntries0 = 5);
  double vals0[] = {2.0, 7.4, 21.0, 8.5, 6.7};
  int inds0[] = {0, 1, 4, 8, 9};
  ECHO(rcpV->ReplaceGlobalValues(NumEntries0, vals0, inds0));

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(int NumEntries = 4);
  double vals[] = {5.1, -2.0, 3.6, 1.1};
  int inds[] = {0, 2, 4, 10};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_SumIntoGlobalValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->SumIntoGlobalValues(NumEntries, vals, inds));

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoMyValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 5);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t mID = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(mID));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(int NumEntries0 = 2);
  double vals0[] = {1.5, 6.9};
  int inds0[] = {1, 4};
  ECHO(rcpV->ReplaceGlobalValues(NumEntries0, vals0, inds0));

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(int NumEntries = 4);
  double vals[] = {-2.0, 5.1};
  int inds[] = {0, 1};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_SumIntoMyValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->SumIntoMyValues(NumEntries, vals, inds));

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

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
double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index );
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

