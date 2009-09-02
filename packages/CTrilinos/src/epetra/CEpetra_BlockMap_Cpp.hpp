#include "CTrilinos_config.h"

#ifndef CEPETRA_BLOCKMAP_CPP_HPP
#define CEPETRA_BLOCKMAP_CPP_HPP


#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const RCP<Epetra_BlockMap>
getBlockMap( CT_Epetra_BlockMap_ID_t id );

/* get Epetra_BlockMap from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_BlockMap>
getBlockMap( CTrilinos_Object_ID_t id );

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const RCP<const Epetra_BlockMap>
getConstBlockMap( CT_Epetra_BlockMap_ID_t id );

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_BlockMap>
getConstBlockMap( CTrilinos_Object_ID_t id );

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
storeBlockMap( Epetra_BlockMap *pobj );

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
storeConstBlockMap( const Epetra_BlockMap *pobj );

/* dump contents of Epetra_BlockMap and const Epetra_BlockMap tables */
void
purgeBlockMapTables(  );


} // namespace CEpetra


#endif // CEPETRA_BLOCKMAP_CPP_HPP


