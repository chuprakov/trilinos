#include "CTrilinos_config.h"

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Export.h"
#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Export>& tableOfExports()
{
    static Table<Epetra_Export>
        loc_tableOfExports(CT_Epetra_Export_ID, "CT_Epetra_Export_ID");
    return loc_tableOfExports;
}


} // namespace


//
// Definitions from CEpetra_Export.h
//


extern "C" {


CT_Epetra_Export_ID_t Epetra_Export_Cast ( 
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        CTrilinos::cast(tableOfExports(), id));
}

CTrilinos_Object_ID_t Epetra_Export_Abstract ( 
  CT_Epetra_Export_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
}

CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pSourceMap = CEpetra::getBlockMap(SourceMapID);
    const Teuchos::RCP<Epetra_BlockMap> 
        pTargetMap = CEpetra::getBlockMap(TargetMapID);

    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store(new Epetra_Export(
        *pSourceMap, *pTargetMap)));

}

CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<Epetra_Export> 
        pExporter = CEpetra::getExport(ExporterID);

    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store(new Epetra_Export(
        *pExporter)));

}

void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID )
{
    CTrilinos_Object_ID_t id =
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(*selfID);
    tableOfExports().remove(&id);
    *selfID = CTrilinos::concreteType<CT_Epetra_Export_ID_t>(id);
}

int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumSameIDs();
}

int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumPermuteIDs();
}

int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->PermuteFromLIDs();
}

int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->PermuteToLIDs();
}

int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumRemoteIDs();
}

int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->RemoteLIDs();
}

int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumExportIDs();
}

int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->ExportLIDs();
}

int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->ExportPIDs();
}

int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumSend();
}

int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getExport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getExport(selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getExport(selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeDistributor(
        &( CEpetra::getExport(selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Export_Cpp.hpp
//


const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CT_Epetra_Export_ID_t id )
{
    return tableOfExports().get(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CTrilinos_Object_ID_t id )
{
    return tableOfExports().get(id);
}

CT_Epetra_Export_ID_t
CEpetra::storeExport( const Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().storeCopy(pobj));
}

void
CEpetra::purgeExportTable(  )
{
    tableOfExports().purge();
}



