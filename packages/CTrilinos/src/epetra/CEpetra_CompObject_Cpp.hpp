#include "CTrilinos_config.h"

#include "CEpetra_CompObject.h"
#include "Epetra_CompObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_COMPOBJECT_CPP_HPP
#define CEPETRA_COMPOBJECT_CPP_HPP


class Epetra_CompObject;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_CompObject from non-const table using CT_Epetra_CompObject_ID */
const RCP<Epetra_CompObject>
getCompObject( CT_Epetra_CompObject_ID_t id );

/* get Epetra_CompObject from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_CompObject>
getCompObject( CTrilinos_Object_ID_t id );

/* get const Epetra_CompObject from either the const or non-const table
 * using CT_Epetra_CompObject_ID */
const RCP<const Epetra_CompObject>
getConstCompObject( CT_Epetra_CompObject_ID_t id );

/* get const Epetra_CompObject from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_CompObject>
getConstCompObject( CTrilinos_Object_ID_t id );

/* store Epetra_CompObject in non-const table */
CT_Epetra_CompObject_ID_t
storeCompObject( Epetra_CompObject *pobj );

/* store const Epetra_CompObject in const table */
CT_Epetra_CompObject_ID_t
storeConstCompObject( const Epetra_CompObject *pobj );

/* dump contents of Epetra_CompObject and const Epetra_CompObject tables */
void
purgeCompObjectTables(  );


} // namespace CEpetra


#endif // CEPETRA_COMPOBJECT_CPP_HPP


