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


/* table to hold objects of type Epetra_Import */
Table<Epetra_Import>& tableOfImports()
{
    static Table<Epetra_Import>
        loc_tableOfImports(CT_Epetra_Import_ID, "CT_Epetra_Import_ID", false);
    return loc_tableOfImports;
}

/* table to hold objects of type const Epetra_Import */
Table<const Epetra_Import>& tableOfConstImports()
{
    static Table<const Epetra_Import>
        loc_tableOfConstImports(CT_Epetra_Import_ID, "CT_Epetra_Import_ID", true);
    return loc_tableOfConstImports;
}


} // namespace


//
// Definitions from CEpetra_Import.h
//


extern "C" {


CT_Epetra_Import_ID_t Epetra_Import_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstImports(), id);
    } else {
        newid = CTrilinos::cast(tableOfImports(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Import_Abstract ( 
  CT_Epetra_Import_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id);
}

CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
            tableOfImports().store(new Epetra_Import(
            *CEpetra::getBlockMap(TargetMapID), *CEpetra::getBlockMap(SourceMapID))));
}

CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
            tableOfImports().store(new Epetra_Import(
            *CEpetra::getImport(ImporterID))));
}

void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Import_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstImports().remove(&aid);
    } else {
        tableOfImports().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Import_ID_t>(aid);
}

int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumSameIDs();
}

int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumPermuteIDs();
}

int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->PermuteFromLIDs();
}

int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->PermuteToLIDs();
}

int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumRemoteIDs();
}

int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->RemoteLIDs();
}

int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumExportIDs();
}

int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->ExportLIDs();
}

int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->ExportPIDs();
}

int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumSend();
}

int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(
        &( CEpetra::getConstImport(selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(
        &( CEpetra::getConstImport(selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeDistributor(
        &( CEpetra::getConstImport(selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Import_Cpp.hpp
//


/* get Epetra_Import from non-const table using CT_Epetra_Import_ID */
const Teuchos::RCP<Epetra_Import>
CEpetra::getImport( CT_Epetra_Import_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id);
    return tableOfImports().get(aid);
}

/* get Epetra_Import from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Import>
CEpetra::getImport( CTrilinos_Object_ID_t id )
{
    return tableOfImports().get(id);
}

/* get const Epetra_Import from either the const or non-const table
 * using CT_Epetra_Import_ID */
const Teuchos::RCP<const Epetra_Import>
CEpetra::getConstImport( CT_Epetra_Import_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id);
    if (id.is_const) {
        return tableOfConstImports().get(aid);
    } else {
        return tableOfImports().get(aid);
    }
}

/* get const Epetra_Import from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Import>
CEpetra::getConstImport( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstImports().get(id);
    } else {
        return tableOfImports().get(id);
    }
}

/* store Epetra_Import in non-const table */
CT_Epetra_Import_ID_t
CEpetra::storeImport( Epetra_Import *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
            tableOfImports().storeCopy(pobj));
}

/* store const Epetra_Import in const table */
CT_Epetra_Import_ID_t
CEpetra::storeConstImport( const Epetra_Import *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
            tableOfConstImports().storeCopy(pobj));
}

/* dump contents of Epetra_Import and const Epetra_Import tables */
void
CEpetra::purgeImportTables(  )
{
    tableOfImports().purge();
    tableOfConstImports().purge();
}



