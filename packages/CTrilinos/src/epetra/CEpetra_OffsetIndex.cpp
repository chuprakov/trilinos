#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_OffsetIndex.h"
#include "Epetra_OffsetIndex.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_OffsetIndex>& tableOfOffsetIndexs()
{
    static Table<Epetra_OffsetIndex> loc_tableOfOffsetIndexs(CT_Epetra_OffsetIndex_ID, "CT_Epetra_OffsetIndex_ID");
    return loc_tableOfOffsetIndexs;
}


} // namespace


//
// Definitions from CEpetra_OffsetIndex.h
//


extern "C" {


CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfOffsetIndexs(), id);
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID )
{
    const Teuchos::RCP<Epetra_CrsGraph> 
        pSourceGraph = CEpetra::getCrsGraph(SourceGraphID);
    const Teuchos::RCP<Epetra_CrsGraph> 
        pTargetGraph = CEpetra::getCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Import> 
        pImporter = CEpetra::getImport(ImporterID);

    return tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *pSourceGraph, *pTargetGraph, *pImporter));

}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<Epetra_CrsGraph> 
        pSourceGraph = CEpetra::getCrsGraph(SourceGraphID);
    const Teuchos::RCP<Epetra_CrsGraph> 
        pTargetGraph = CEpetra::getCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Export> 
        pExporter = CEpetra::getExport(ExporterID);

    return tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *pSourceGraph, *pTargetGraph, *pExporter));

}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<Epetra_OffsetIndex> 
        pIndexor = CEpetra::getOffsetIndex(IndexorID);

    return tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *pIndexor));

}

void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID )
{
    tableOfOffsetIndexs().remove(selfID);
}

int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getOffsetIndex(selfID)->SameOffsets();
}

int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getOffsetIndex(selfID)->PermuteOffsets();
}

int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getOffsetIndex(selfID)->RemoteOffsets();
}


} // extern "C"


//
// Definitions from CEpetra_OffsetIndex_Cpp.hpp
//


const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    return tableOfOffsetIndexs().get(id);
}

CT_Epetra_OffsetIndex_ID_t
CEpetra::storeOffsetIndex( const Epetra_OffsetIndex *pobj )
{
    return tableOfOffsetIndexs().store(pobj);
}

void
CEpetra::purgeOffsetIndexTable(  )
{
    tableOfOffsetIndexs().purge();
}



