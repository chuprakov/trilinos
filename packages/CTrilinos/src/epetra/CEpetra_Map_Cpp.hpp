#include "CTrilinos_config.h"

#ifndef CEPETRA_MAP_CPP_HPP
#define CEPETRA_MAP_CPP_HPP


#include "CEpetra_Map.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_Map;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const RCP<Epetra_Map>
getMap( CT_Epetra_Map_ID_t id );

/* get Epetra_Map from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Map>
getMap( CTrilinos_Object_ID_t id );

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const RCP<const Epetra_Map>
getConstMap( CT_Epetra_Map_ID_t id );

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Map>
getConstMap( CTrilinos_Object_ID_t id );

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
storeMap( Epetra_Map *pobj );

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
storeConstMap( const Epetra_Map *pobj );

/* dump contents of Epetra_Map and const Epetra_Map tables */
void
purgeMapTables(  );


} // namespace CEpetra


#endif // CEPETRA_MAP_CPP_HPP


