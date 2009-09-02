#include "CTrilinos_config.h"

#ifndef CEPETRA_SRCDISTOBJECT_CPP_HPP
#define CEPETRA_SRCDISTOBJECT_CPP_HPP


#include "Epetra_SrcDistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_SrcDistObject from non-const table using CT_Epetra_SrcDistObject_ID */
const RCP<Epetra_SrcDistObject>
getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id );

/* get Epetra_SrcDistObject from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_SrcDistObject>
getSrcDistObject( CTrilinos_Object_ID_t id );

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CT_Epetra_SrcDistObject_ID */
const RCP<const Epetra_SrcDistObject>
getConstSrcDistObject( CT_Epetra_SrcDistObject_ID_t id );

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_SrcDistObject>
getConstSrcDistObject( CTrilinos_Object_ID_t id );

/* store Epetra_SrcDistObject in non-const table */
CT_Epetra_SrcDistObject_ID_t
storeSrcDistObject( Epetra_SrcDistObject *pobj );

/* store const Epetra_SrcDistObject in const table */
CT_Epetra_SrcDistObject_ID_t
storeConstSrcDistObject( const Epetra_SrcDistObject *pobj );

/* dump contents of Epetra_SrcDistObject and const Epetra_SrcDistObject tables */
void
purgeSrcDistObjectTables(  );


} // namespace CEpetra


#endif // CEPETRA_SRCDISTOBJECT_CPP_HPP


