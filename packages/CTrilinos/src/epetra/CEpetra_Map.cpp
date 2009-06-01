#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Map.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Map>& tableOfMaps()
{
    static Table<Epetra_Map> loc_tableOfMaps(CT_Epetra_Map_ID, "CT_Epetra_Map_ID");
    return loc_tableOfMaps;
}


} // namespace


//
// Definitions from CEpetra_Map.h
//


extern "C" {


CT_Epetra_Map_ID_t Epetra_Map_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfMaps(), id);
}

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<Epetra_Comm> 
        pComm = CEpetra::getComm(CommID);

    return tableOfMaps().store(new Epetra_Map(
        NumGlobalElements, IndexBase, *pComm));

}

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<Epetra_Comm> 
        pComm = CEpetra::getComm(CommID);

    return tableOfMaps().store(new Epetra_Map(
        NumGlobalElements, NumMyElements, IndexBase, *pComm));

}

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<Epetra_Comm> 
        pComm = CEpetra::getComm(CommID);

    return tableOfMaps().store(new Epetra_Map(
        NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, *pComm));

}

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID )
{
    const Teuchos::RCP<Epetra_Map> 
        pmap = CEpetra::getMap(mapID);

    return tableOfMaps().store(new Epetra_Map(*pmap));

}

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID )
{
    tableOfMaps().remove(selfID);
}


} // extern "C"


//
// Definitions from CEpetra_Map_Cpp.hpp
//


const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CT_Epetra_Map_ID_t id )
{
    return tableOfMaps().get(id);
}

void
CEpetra::purgeMapTable(  )
{
    tableOfMaps().purge();
}



