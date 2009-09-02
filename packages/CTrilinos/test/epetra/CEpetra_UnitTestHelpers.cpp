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
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CAmesos_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CTeuchos_any_Cpp.hpp"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"

#include "CEpetra_SerialComm.h"
#include "CEpetra_Comm.h"

#include "CEpetra_UnitTestHelpers.hpp"

void CEpetra_Test_CleanSlate()
{
  CEpetra::purgeCommTables();
  CEpetra::purgeSerialCommTables();
#ifdef HAVE_MPI
  CEpetra::purgeMpiCommTables();
#endif /* HAVE_MPI */
  CEpetra::purgeBlockMapTables();
  CEpetra::purgeMapTables();
  CEpetra::purgeObjectTables();
  CEpetra::purgeDistributorTables();
  CEpetra::purgeDirectoryTables();
  CEpetra::purgeImportTables();
  CEpetra::purgeExportTables();
  CEpetra::purgeOffsetIndexTables();
  CEpetra::purgeDistObjectTables();
  CEpetra::purgeSrcDistObjectTables();
  CEpetra::purgeBLASTables();
  CEpetra::purgeFlopsTables();
  CEpetra::purgeCompObjectTables();
  CEpetra::purgeMultiVectorTables();
  CEpetra::purgeVectorTables();
  CEpetra::purgeCrsGraphTables();
  CEpetra::purgeCrsMatrixTables();
  CEpetra::purgeOperatorTables();
  CEpetra::purgeRowMatrixTables();
  CEpetra::purgeTimeTables();
  CEpetra::purgeJadMatrixTables();
  CEpetra::purgeLinearProblemTables();
  CEpetra::purgeLAPACKTables();
  CAmesos::purgeAmesosTables();
  CAmesos::purgeBaseSolverTables();
  CTeuchos::purgeanyTables();
  CTeuchos::purgeCommandLineProcessorTables();
  CTeuchos::purgeParameterEntryTables();
  CTeuchos::purgeParameterListTables();
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

