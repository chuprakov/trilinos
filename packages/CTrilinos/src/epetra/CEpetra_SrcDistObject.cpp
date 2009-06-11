#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_SrcDistObject.h"
#include "Epetra_SrcDistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_SrcDistObject>& tableOfSrcDistObjects()
{
    static Table<Epetra_SrcDistObject> loc_tableOfSrcDistObjects(CT_Epetra_SrcDistObject_ID, "CT_Epetra_SrcDistObject_ID");
    return loc_tableOfSrcDistObjects;
}


} // namespace


//
// Definitions from CEpetra_SrcDistObject.h
//


extern "C" {


CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfSrcDistObjects(), id);
}

void Epetra_SrcDistObject_Destroy ( 
  CT_Epetra_SrcDistObject_ID_t * selfID )
{
    tableOfSrcDistObjects().remove(selfID);
}

CT_Epetra_BlockMap_ID_t Epetra_SrcDistObject_Map ( 
  CT_Epetra_SrcDistObject_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getSrcDistObject(selfID)->Map() ));
}


} // extern "C"


//
// Definitions from CEpetra_SrcDistObject_Cpp.hpp
//


const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    return tableOfSrcDistObjects().get(id);
}

CT_Epetra_SrcDistObject_ID_t
CEpetra::storeSrcDistObject( const Epetra_SrcDistObject *pobj )
{
    return tableOfSrcDistObjects().store(pobj);
}

void
CEpetra::purgeSrcDistObjectTable(  )
{
    tableOfSrcDistObjects().purge();
}



