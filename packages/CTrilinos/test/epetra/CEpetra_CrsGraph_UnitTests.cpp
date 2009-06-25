#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Export.h"
#include "CEpetra_Import.h"
#include "Epetra_DataAccess.h"
#include "Epetra_CrsGraph.h"
#include "CEpetra_CrsGraph.h"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID );
 **********************************************************************/

/**********************************************************************
void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_NoDiagonal ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_HaveColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_MaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );
 **********************************************************************/

/**********************************************************************
int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc );
 **********************************************************************/

/**********************************************************************
void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );
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

