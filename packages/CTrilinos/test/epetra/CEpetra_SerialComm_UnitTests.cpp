#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Directory.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_SerialCommData.h"
#include "Epetra_SerialComm.h"
#include "CEpetra_SerialComm.h"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* This cast should be allowed */
  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());
  ECHO(CT_Epetra_SerialComm_ID_t dupID = Epetra_SerialComm_Cast(
       selfID));

  /* This cast should be allowed */
  ECHO(CT_Epetra_Comm_ID_t commID = Epetra_Comm_Cast(selfID));
  ECHO(CT_Epetra_SerialComm_ID_t serialcommID = Epetra_SerialComm_Cast(
       commID));

  /* If no exceptions thrown, then test was successful */
}

/**********************************************************************
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_SerialComm_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(CT_Epetra_SerialComm_ID_t dupID = Epetra_SerialComm_Duplicate(
       selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_SerialComm_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
}

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_SerialComm_Destroy ( 
  CT_Epetra_SerialComm_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(Epetra_SerialComm_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MyPID )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(int ret = Epetra_SerialComm_MyPID(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , NumProc )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(int ret = Epetra_SerialComm_NumProc(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(ret, 1);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_SerialCommData_ID_t Epetra_SerialComm_DataPtr ( 
  CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_SerialComm_ID_t CommID );
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

