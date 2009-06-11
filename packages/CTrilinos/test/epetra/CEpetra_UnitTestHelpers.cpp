#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm_Cpp.hpp"
#endif /* HAVE_MPI */
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"

#include "CEpetra_UnitTestHelpers.hpp"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )


#define BASE_COMM Epetra_Comm
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


void CEpetra_Test_CleanSlate()
{
  CEpetra::purgeCommTable();
  CEpetra::purgeMultiVectorTable();
  CEpetra::purgeMapTable();
#ifdef HAVE_MPI
  CEpetra::purgeMpiCommTable();
#endif /* HAVE_MPI */
  CEpetra::purgeBlockMapTable();
  CEpetra::purgeSerialCommTable();
  CEpetra::purgeVectorTable();
}


CT_Epetra_Comm_ID_t
UnitTest_Create_Comm()
{
  return Epetra_Comm_Cast(COMM_CREATE ());
}

CT_Epetra_Map_ID_t
UnitTest_Create_Map(int NumGlobalElements)
{
  return Epetra_Map_Create(NumGlobalElements, 0, UnitTest_Create_Comm());
}

CT_Epetra_BlockMap_ID_t
UnitTest_Create_BlockMapCast(int NumGlobalElements)
{
  return Epetra_BlockMap_Cast(UnitTest_Create_Map(NumGlobalElements));
}

CT_Epetra_Vector_ID_t
UnitTest_Create_Vector(int NumGlobalElements)
{
  return Epetra_Vector_Create(UnitTest_Create_BlockMapCast(NumGlobalElements), true);
}

