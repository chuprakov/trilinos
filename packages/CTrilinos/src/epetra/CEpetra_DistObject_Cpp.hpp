#include "CTrilinos_config.h"

#ifndef CEPETRA_DISTOBJECT_CPP_HPP
#define CEPETRA_DISTOBJECT_CPP_HPP


#include "Epetra_DistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const RCP<Epetra_DistObject>
getDistObject( CT_Epetra_DistObject_ID_t id );

/* get Epetra_DistObject from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_DistObject>
getDistObject( CTrilinos_Object_ID_t id );

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const RCP<const Epetra_DistObject>
getConstDistObject( CT_Epetra_DistObject_ID_t id );

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_DistObject>
getConstDistObject( CTrilinos_Object_ID_t id );

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
storeDistObject( Epetra_DistObject *pobj );

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
storeConstDistObject( const Epetra_DistObject *pobj );

/* dump contents of Epetra_DistObject and const Epetra_DistObject tables */
void
purgeDistObjectTables(  );


} // namespace CEpetra


#endif // CEPETRA_DISTOBJECT_CPP_HPP


