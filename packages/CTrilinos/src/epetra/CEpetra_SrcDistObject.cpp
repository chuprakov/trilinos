#include "CTrilinos_config.h"

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_SrcDistObject.h"
#include "Epetra_SrcDistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_SrcDistObject */
Table<Epetra_SrcDistObject>& tableOfSrcDistObjects()
{
    static Table<Epetra_SrcDistObject>
        loc_tableOfSrcDistObjects(CT_Epetra_SrcDistObject_ID, "CT_Epetra_SrcDistObject_ID", false);
    return loc_tableOfSrcDistObjects;
}

/* table to hold objects of type const Epetra_SrcDistObject */
Table<const Epetra_SrcDistObject>& tableOfConstSrcDistObjects()
{
    static Table<const Epetra_SrcDistObject>
        loc_tableOfConstSrcDistObjects(CT_Epetra_SrcDistObject_ID, "CT_Epetra_SrcDistObject_ID", true);
    return loc_tableOfConstSrcDistObjects;
}


} // namespace


//
// Definitions from CEpetra_SrcDistObject.h
//


extern "C" {


CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstSrcDistObjects(), id);
    } else {
        newid = CTrilinos::cast(tableOfSrcDistObjects(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_SrcDistObject_Abstract ( 
  CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id);
}

void Epetra_SrcDistObject_Destroy ( 
  CT_Epetra_SrcDistObject_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstSrcDistObjects().remove(&aid);
    } else {
        tableOfSrcDistObjects().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(aid);
}

CT_Epetra_BlockMap_ID_t Epetra_SrcDistObject_Map ( 
  CT_Epetra_SrcDistObject_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(
        &( CEpetra::getConstSrcDistObject(selfID)->Map() ));
}


} // extern "C"


//
// Definitions from CEpetra_SrcDistObject_Cpp.hpp
//


/* get Epetra_SrcDistObject from non-const table using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id);
    return tableOfSrcDistObjects().get(aid);
}

/* get Epetra_SrcDistObject from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CTrilinos_Object_ID_t id )
{
    return tableOfSrcDistObjects().get(id);
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id);
    if (id.is_const) {
        return tableOfConstSrcDistObjects().get(aid);
    } else {
        return tableOfSrcDistObjects().get(aid);
    }
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstSrcDistObjects().get(id);
    } else {
        return tableOfSrcDistObjects().get(id);
    }
}

/* store Epetra_SrcDistObject in non-const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeSrcDistObject( Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
            tableOfSrcDistObjects().storeCopy(pobj));
}

/* store const Epetra_SrcDistObject in const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeConstSrcDistObject( const Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
            tableOfConstSrcDistObjects().storeCopy(pobj));
}

/* dump contents of Epetra_SrcDistObject and const Epetra_SrcDistObject tables */
void
CEpetra::purgeSrcDistObjectTables(  )
{
    tableOfSrcDistObjects().purge();
    tableOfConstSrcDistObjects().purge();
}



