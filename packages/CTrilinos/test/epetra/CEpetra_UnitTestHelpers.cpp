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
