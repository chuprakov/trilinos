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


Table<Epetra_Directory>& tableOfDirectorys()
{
    static Table<Epetra_Directory> loc_tableOfDirectorys(CT_Epetra_Directory_ID, "CT_Epetra_Directory_ID");
    return loc_tableOfDirectorys;
}


} // namespace


//
// Definitions from CEpetra_Directory.h
//


extern "C" {


CT_Epetra_Directory_ID_t Epetra_Directory_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfDirectorys(), id);
}

void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID )
{
    tableOfDirectorys().remove(selfID);
}

int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return CEpetra::getDirectory(selfID)->GetDirectoryEntries(
        *pMap, NumEntries, GlobalEntries, Procs, LocalEntries, EntrySizes, high_rank_sharing_procs);
}

boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID )
{
    return CEpetra::getDirectory(selfID)->GIDsAllUniquelyOwned();
}


} // extern "C"


//
// Definitions from CEpetra_Directory_Cpp.hpp
//


const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CT_Epetra_Directory_ID_t id )
{
    return tableOfDirectorys().get(id);
}

CT_Epetra_Directory_ID_t
CEpetra::storeDirectory( const Epetra_Directory *pobj )
{
    return tableOfDirectorys().store(pobj);
}

void
CEpetra::purgeDirectoryTable(  )
{
    tableOfDirectorys().purge();
}



