#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm_Cpp.hpp"
#endif /* HAVE_MPI */
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"

#include "CEpetra_UnitTestHelpers.hpp"

void CEpetra_Test_CleanSlate()
{
  CEpetra::purgeCommTable();
  CEpetra::purgeSerialCommTable();
#ifdef HAVE_MPI
  CEpetra::purgeMpiCommTable();
#endif /* HAVE_MPI */
  CEpetra::purgeBlockMapTable();
  CEpetra::purgeMapTable();
  CEpetra::purgeObjectTable();
  CEpetra::purgeDistributorTable();
  CEpetra::purgeDirectoryTable();
  CEpetra::purgeImportTable();
  CEpetra::purgeExportTable();
  CEpetra::purgeOffsetIndexTable();
  CEpetra::purgeDistObjectTable();
  CEpetra::purgeSrcDistObjectTable();
  CEpetra::purgeBLASTable();
  CEpetra::purgeFlopsTable();
  CEpetra::purgeCompObjectTable();
  CEpetra::purgeMultiVectorTable();
  CEpetra::purgeVectorTable();
  CEpetra::purgeCrsGraphTable();
  CEpetra::purgeCrsMatrixTable();
  CEpetra::purgeOperatorTable();
  CEpetra::purgeRowMatrixTable();
}

CT_Epetra_Comm_ID_t
UnitTest_Create_Comm()
{
#ifdef EPETRA_MPI
  return Epetra_Comm_Cast(Epetra_MpiComm_Abstract(Epetra_MpiComm_Create(MPI_COMM_WORLD)));
#else /* EPETRA_MPI */
  return Epetra_Comm_Cast(Epetra_SerialComm_Abstract(Epetra_SerialComm_Create()));
#endif /* EPETRA_MPI */
}

