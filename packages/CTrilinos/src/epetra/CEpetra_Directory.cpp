#include "CTrilinos_config.h"

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Directory.h"
#include "Epetra_Directory.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Directory */
Table<Epetra_Directory>& tableOfDirectorys()
{
    static Table<Epetra_Directory>
        loc_tableOfDirectorys(CT_Epetra_Directory_ID, "CT_Epetra_Directory_ID", false);
    return loc_tableOfDirectorys;
}

/* table to hold objects of type const Epetra_Directory */
Table<const Epetra_Directory>& tableOfConstDirectorys()
{
    static Table<const Epetra_Directory>
        loc_tableOfConstDirectorys(CT_Epetra_Directory_ID, "CT_Epetra_Directory_ID", true);
    return loc_tableOfConstDirectorys;
}


} // namespace


//
// Definitions from CEpetra_Directory.h
//


extern "C" {


CT_Epetra_Directory_ID_t Epetra_Directory_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstDirectorys(), id);
    } else {
        newid = CTrilinos::cast(tableOfDirectorys(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Directory_Abstract ( 
  CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id);
}

void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstDirectorys().remove(&aid);
    } else {
        tableOfDirectorys().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(aid);
}

int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs )
{
    return CEpetra::getConstDirectory(selfID)->GetDirectoryEntries(
            *CEpetra::getBlockMap(MapID), NumEntries, GlobalEntries, Procs, LocalEntries, EntrySizes, high_rank_sharing_procs);
}

boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID )
{
    return CEpetra::getConstDirectory(selfID)->GIDsAllUniquelyOwned();
}


} // extern "C"


//
// Definitions from CEpetra_Directory_Cpp.hpp
//


/* get Epetra_Directory from non-const table using CT_Epetra_Directory_ID */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CT_Epetra_Directory_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id);
    return tableOfDirectorys().get(aid);
}

/* get Epetra_Directory from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CTrilinos_Object_ID_t id )
{
    return tableOfDirectorys().get(id);
}

/* get const Epetra_Directory from either the const or non-const table
 * using CT_Epetra_Directory_ID */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CT_Epetra_Directory_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id);
    if (id.is_const) {
        return tableOfConstDirectorys().get(aid);
    } else {
        return tableOfDirectorys().get(aid);
    }
}

/* get const Epetra_Directory from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstDirectorys().get(id);
    } else {
        return tableOfDirectorys().get(id);
    }
}

/* store Epetra_Directory in non-const table */
CT_Epetra_Directory_ID_t
CEpetra::storeDirectory( Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
            tableOfDirectorys().storeCopy(pobj));
}

/* store const Epetra_Directory in const table */
CT_Epetra_Directory_ID_t
CEpetra::storeConstDirectory( const Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
            tableOfConstDirectorys().storeCopy(pobj));
}

/* dump contents of Epetra_Directory and const Epetra_Directory tables */
void
CEpetra::purgeDirectoryTables(  )
{
    tableOfDirectorys().purge();
    tableOfConstDirectorys().purge();
}



