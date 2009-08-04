#include "CTrilinos_config.h"

#ifndef CEPETRA_OFFSETINDEX_CPP_HPP
#define CEPETRA_OFFSETINDEX_CPP_HPP


#include "CEpetra_OffsetIndex.h"
#include "Epetra_OffsetIndex.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_OffsetIndex;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_OffsetIndex from non-const table using CT_Epetra_OffsetIndex_ID */
const RCP<Epetra_OffsetIndex>
getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id );

/* get Epetra_OffsetIndex from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_OffsetIndex>
getOffsetIndex( CTrilinos_Object_ID_t id );

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CT_Epetra_OffsetIndex_ID */
const RCP<const Epetra_OffsetIndex>
getConstOffsetIndex( CT_Epetra_OffsetIndex_ID_t id );

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_OffsetIndex>
getConstOffsetIndex( CTrilinos_Object_ID_t id );

/* store Epetra_OffsetIndex in non-const table */
CT_Epetra_OffsetIndex_ID_t
storeOffsetIndex( Epetra_OffsetIndex *pobj );

/* store const Epetra_OffsetIndex in const table */
CT_Epetra_OffsetIndex_ID_t
storeConstOffsetIndex( const Epetra_OffsetIndex *pobj );

/* dump contents of Epetra_OffsetIndex and const Epetra_OffsetIndex tables */
void
purgeOffsetIndexTables(  );


} // namespace CEpetra


#endif // CEPETRA_OFFSETINDEX_CPP_HPP


