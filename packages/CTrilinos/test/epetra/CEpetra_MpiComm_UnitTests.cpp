#ifdef HAVE_MPI

#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Directory.h"
#include "CEpetra_Distributor.h"
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "CEpetra_MpiComm.h"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* This cast should be allowed */
  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));
  ECHO(CT_Epetra_MpiComm_ID_t dupID = Epetra_MpiComm_Cast(Epetra_MpiComm_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), true);

  /* These casts should be allowed */
  ECHO(CT_Epetra_Comm_ID_t commID = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, commID), true);
  ECHO(CT_Epetra_MpiComm_ID_t mpicommID = Epetra_MpiComm_Cast(Epetra_Comm_Abstract(commID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(commID, mpicommID), true);

  /* If no exceptions thrown, then test was successful */
}

/**********************************************************************
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_MpiComm_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_MpiComm_ID_t dupID = Epetra_MpiComm_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_MpiComm_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Clone )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_Comm_ID_t dupID = Epetra_MpiComm_Clone(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_Comm_ID);
  TEST_EQUALITY_CONST(dupID.index, 0);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Barrier )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Barrier(selfID));

  /* ??? */
}

/**********************************************************************
int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Comm )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(MPI_Comm c = Epetra_MpiComm_Comm(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */
}

/**********************************************************************
int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , MyPID )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(int ret = Epetra_MpiComm_MyPID(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */
}

/**********************************************************************
int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , NumProc )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(int ret = Epetra_MpiComm_NumProc(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , CreateDistributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_Distributor_ID_t disID = Epetra_MpiComm_CreateDistributor(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(disID.type, CT_Epetra_Distributor_ID);
  TEST_EQUALITY_CONST(disID.index, 0);
}

/**********************************************************************
CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , CreateDirectory )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));
  ECHO(CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Cast(Epetra_MpiComm_Abstract(selfID)));

  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(CT_Epetra_Directory_ID_t dirID = Epetra_MpiComm_CreateDirectory(selfID, MapID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dirID.type, CT_Epetra_Directory_ID);
  TEST_EQUALITY_CONST(dirID.index, 0);
}

/**********************************************************************
int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t commID = Epetra_MpiComm_Create(MPI_COMM_WORLD));
  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Assign(selfID, commID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, commID), false);
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

#endif /* HAVE_MPI */

