#include "CTrilinos_config.h"

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Import.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Import>& tableOfImports()
{
    static Table<Epetra_Import> loc_tableOfImports(CT_Epetra_Import_ID, "CT_Epetra_Import_ID");
    return loc_tableOfImports;
}


} // namespace


//
// Definitions from CEpetra_Import.h
//


extern "C" {


CT_Epetra_Import_ID_t Epetra_Import_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfImports(), id);
}

CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pTargetMap = CEpetra::getBlockMap(TargetMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pSourceMap = CEpetra::getBlockMap(SourceMapID);

    return tableOfImports().store(new Epetra_Import(
        *pTargetMap, *pSourceMap));

}

CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID )
{
    const Teuchos::RCP<Epetra_Import> 
        pImporter = CEpetra::getImport(ImporterID);

    return tableOfImports().store(new Epetra_Import(*pImporter));

}

void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID )
{
    tableOfImports().remove(selfID);
}

int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumSameIDs();
}

int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumPermuteIDs();
}

int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->PermuteFromLIDs();
}

int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->PermuteToLIDs();
}

int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumRemoteIDs();
}

int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->RemoteLIDs();
}

int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumExportIDs();
}

int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->ExportLIDs();
}

int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->ExportPIDs();
}

int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumSend();
}

int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getImport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getImport(selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getImport(selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeDistributor(
        &( CEpetra::getImport(selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Import_Cpp.hpp
//


const Teuchos::RCP<Epetra_Import>
CEpetra::getImport( CT_Epetra_Import_ID_t id )
{
    return tableOfImports().get(id);
}

CT_Epetra_Import_ID_t
CEpetra::storeImport( const Epetra_Import *pobj )
{
    return tableOfImports().storeCopy(pobj);
}

void
CEpetra::purgeImportTable(  )
{
    tableOfImports().purge();
}



