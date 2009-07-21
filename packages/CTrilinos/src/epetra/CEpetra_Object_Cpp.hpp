#include "CTrilinos_config.h"

#include "CEpetra_Object.h"
#include "Epetra_Object.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_OBJECT_CPP_HPP
#define CEPETRA_OBJECT_CPP_HPP


class Epetra_Object;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Object from non-const table using CT_Epetra_Object_ID */
const RCP<Epetra_Object>
getObject( CT_Epetra_Object_ID_t id );

/* get Epetra_Object from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Object>
getObject( CTrilinos_Object_ID_t id );

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const RCP<const Epetra_Object>
getConstObject( CT_Epetra_Object_ID_t id );

/* get const Epetra_Object from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Object>
getConstObject( CTrilinos_Object_ID_t id );

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
storeObject( Epetra_Object *pobj );

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
storeConstObject( const Epetra_Object *pobj );

/* dump contents of Epetra_Object and const Epetra_Object tables */
void
purgeObjectTables(  );


} // namespace CEpetra


#endif // CEPETRA_OBJECT_CPP_HPP


