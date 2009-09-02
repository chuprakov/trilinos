#include "CTrilinos_config.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Map.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Map */
Table<Epetra_Map>& tableOfMaps()
{
    static Table<Epetra_Map>
        loc_tableOfMaps(CT_Epetra_Map_ID, "CT_Epetra_Map_ID", false);
    return loc_tableOfMaps;
}

/* table to hold objects of type const Epetra_Map */
Table<const Epetra_Map>& tableOfConstMaps()
{
    static Table<const Epetra_Map>
        loc_tableOfConstMaps(CT_Epetra_Map_ID, "CT_Epetra_Map_ID", true);
    return loc_tableOfConstMaps;
}


} // namespace


//
// Definitions from CEpetra_Map.h
//


extern "C" {


CT_Epetra_Map_ID_t Epetra_Map_Cast ( CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstMaps(), id);
    } else {
        newid = CTrilinos::cast(tableOfMaps(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Map_Abstract ( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
}

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store(new Epetra_Map(NumGlobalElements, 
        IndexBase, *CEpetra::getConstComm(CommID))));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store(new Epetra_Map(NumGlobalElements, 
        NumMyElements, IndexBase, *CEpetra::getConstComm(CommID))));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store(new Epetra_Map(NumGlobalElements, 
        NumMyElements, MyGlobalElements, IndexBase, 
        *CEpetra::getConstComm(CommID))));
}

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store(new Epetra_Map(*CEpetra::getConstMap(
        mapID))));
}

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Map_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstMaps().remove(&aid);
    } else {
        tableOfMaps().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Map_ID_t>(aid);
}

void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID )
{
    Epetra_Map& self = *( CEpetra::getMap(selfID) );

    self = *CEpetra::getConstMap(mapID);
}


} // extern "C"


//
// Definitions from CEpetra_Map_Cpp.hpp
//


/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CT_Epetra_Map_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
    return tableOfMaps().get(aid);
}

/* get Epetra_Map from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CTrilinos_Object_ID_t id )
{
    return tableOfMaps().get(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CT_Epetra_Map_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
    if (id.is_const) {
        return tableOfConstMaps().get(aid);
    } else {
        return tableOfMaps().get(aid);
    }
}

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstMaps().get(id);
    } else {
        return tableOfMaps().get(id);
    }
}

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeMap( Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
            tableOfMaps().storeShared(pobj));
}

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
CEpetra::storeConstMap( const Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
            tableOfConstMaps().storeShared(pobj));
}

/* dump contents of Epetra_Map and const Epetra_Map tables */
void
CEpetra::purgeMapTables(  )
{
    tableOfMaps().purge();
    tableOfConstMaps().purge();
}



