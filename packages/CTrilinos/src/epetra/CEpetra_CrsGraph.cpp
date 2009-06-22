#include "CEpetra_BlockMap_Cpp.hpp"
#include "Epetra_DataAccess.h"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_CrsGraphData_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_CrsGraph.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_CrsGraph>& tableOfCrsGraphs()
{
    static Table<Epetra_CrsGraph> loc_tableOfCrsGraphs(CT_Epetra_CrsGraph_ID, "CT_Epetra_CrsGraph_ID");
    return loc_tableOfCrsGraphs;
}


} // namespace


//
// Definitions from CEpetra_CrsGraph.h
//


extern "C" {


CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfCrsGraphs(), id);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pRowMap = CEpetra::getBlockMap(RowMapID);

    return tableOfCrsGraphs().store(new Epetra_CrsGraph(
        CV, *pRowMap, NumIndicesPerRow, StaticProfile));

}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pRowMap = CEpetra::getBlockMap(RowMapID);

    return tableOfCrsGraphs().store(new Epetra_CrsGraph(
        CV, *pRowMap, NumIndicesPerRow, StaticProfile));

}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pRowMap = CEpetra::getBlockMap(RowMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pColMap = CEpetra::getBlockMap(ColMapID);

    return tableOfCrsGraphs().store(new Epetra_CrsGraph(
        CV, *pRowMap, *pColMap, NumIndicesPerRow, StaticProfile));

}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pRowMap = CEpetra::getBlockMap(RowMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pColMap = CEpetra::getBlockMap(ColMapID);

    return tableOfCrsGraphs().store(new Epetra_CrsGraph(
        CV, *pRowMap, *pColMap, NumIndicesPerRow, StaticProfile));

}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<Epetra_CrsGraph> 
        pGraph = CEpetra::getCrsGraph(GraphID);

    return tableOfCrsGraphs().store(new Epetra_CrsGraph(*pGraph));

}

void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID )
{
    tableOfCrsGraphs().remove(selfID);
}

int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertGlobalIndices(
        GlobalRow, NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(
        GlobalRow, NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(Row);
}

int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertMyIndices(
        LocalRow, NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(
        LocalRow, NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(Row);
}

int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->FillComplete();
}

int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pDomainMap = CEpetra::getBlockMap(DomainMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pRangeMap = CEpetra::getBlockMap(RangeMapID);

    return CEpetra::getCrsGraph(selfID)->FillComplete(
        *pDomainMap, *pRangeMap);
}

int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->OptimizeStorage();
}

int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->ExtractGlobalRowCopy(
        GlobalRow, LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->ExtractMyRowCopy(
        LocalRow, LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getCrsGraph(selfID)->ExtractGlobalRowView(
        GlobalRow, *NumIndices, *Indices);
}

int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getCrsGraph(selfID)->ExtractMyRowView(
        LocalRow, *NumIndices, *Indices);
}

boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->Filled();
}

boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->StorageOptimized();
}

boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->IndicesAreGlobal();
}

boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->IndicesAreLocal();
}

boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->LowerTriangular();
}

boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->UpperTriangular();
}

boolean Epetra_CrsGraph_NoDiagonal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NoDiagonal();
}

boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID )
{
    return CEpetra::getCrsGraph(selfID)->MyGlobalRow(GID);
}

boolean Epetra_CrsGraph_HaveColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->HaveColMap();
}

int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyRows();
}

int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalRows();
}

int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyCols();
}

int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalCols();
}

int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsGraph_NumMyDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyDiagonals();
}

int Epetra_CrsGraph_NumMyBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyBlockRows();
}

int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalBlockRows();
}

int Epetra_CrsGraph_NumMyBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyBlockCols();
}

int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalBlockCols();
}

int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalEntries();
}

int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyEntries();
}

int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->MaxRowDim();
}

int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->GlobalMaxRowDim();
}

int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->MaxColDim();
}

int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->GlobalMaxColDim();
}

int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->NumMyNonzeros();
}

int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->NumGlobalIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->NumAllocatedGlobalIndices(
        Row);
}

int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->MaxNumIndices();
}

int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->GlobalMaxNumIndices();
}

int Epetra_CrsGraph_MaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->MaxNumNonzeros();
}

int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->GlobalMaxNumNonzeros();
}

int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->NumMyIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->NumAllocatedMyIndices(Row);
}

int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->IndexBase();
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getCrsGraph(selfID)->RowMap() ));
}

int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pnewmap = CEpetra::getBlockMap(newmapID);

    return CEpetra::getCrsGraph(selfID)->ReplaceRowMap(*pnewmap);
}

int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pnewmap = CEpetra::getBlockMap(newmapID);

    return CEpetra::getCrsGraph(selfID)->ReplaceColMap(*pnewmap);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getCrsGraph(selfID)->ColMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getCrsGraph(selfID)->DomainMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getCrsGraph(selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeImport(
        CEpetra::getCrsGraph(selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeExport(
        CEpetra::getCrsGraph(selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeComm(
        &( CEpetra::getCrsGraph(selfID)->Comm() ));
}

int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return CEpetra::getCrsGraph(selfID)->LRID(GRID_in);
}

int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return CEpetra::getCrsGraph(selfID)->GRID(LRID_in);
}

int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return CEpetra::getCrsGraph(selfID)->LCID(GCID_in);
}

int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return CEpetra::getCrsGraph(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return CEpetra::getCrsGraph(selfID)->MyGRID(GRID_in);
}

boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return CEpetra::getCrsGraph(selfID)->MyLRID(LRID_in);
}

boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return CEpetra::getCrsGraph(selfID)->MyGCID(GCID_in);
}

boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return CEpetra::getCrsGraph(selfID)->MyLCID(LCID_in);
}

int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc )
{
    Epetra_CrsGraph& self = *( CEpetra::getCrsGraph(selfID) );

    return self[Loc];
}

void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID )
{
    Epetra_CrsGraph& self = *( CEpetra::getCrsGraph(selfID) );

    const Teuchos::RCP<Epetra_CrsGraph> 
        pSource = CEpetra::getCrsGraph(SourceID);

    self = *pSource;
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getCrsGraph(selfID)->ImportMap() ));
}

int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->TransformToLocal();
}

int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pDomainMap = CEpetra::getBlockMap(DomainMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pRangeMap = CEpetra::getBlockMap(RangeMapID);

    return CEpetra::getCrsGraph(selfID)->TransformToLocal(
        &(*pDomainMap), &(*pRangeMap));
}

CT_Epetra_CrsGraphData_ID_t Epetra_CrsGraph_DataPtr ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeCrsGraphData(
        CEpetra::getCrsGraph(selfID)->DataPtr());
}

void Epetra_CrsGraph_SortGhostsAssociatedWithEachProcessor ( 
  CT_Epetra_CrsGraph_ID_t selfID, boolean Flag )
{
    CEpetra::getCrsGraph(selfID)->SortGhostsAssociatedWithEachProcessor(
        Flag);
}


} // extern "C"


//
// Definitions from CEpetra_CrsGraph_Cpp.hpp
//


const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return tableOfCrsGraphs().get(id);
}

CT_Epetra_CrsGraph_ID_t
CEpetra::storeCrsGraph( const Epetra_CrsGraph *pobj )
{
    return tableOfCrsGraphs().store(pobj);
}

void
CEpetra::purgeCrsGraphTable(  )
{
    tableOfCrsGraphs().purge();
}



