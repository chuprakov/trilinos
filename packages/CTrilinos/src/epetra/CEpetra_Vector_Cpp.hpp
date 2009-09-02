#include "CTrilinos_config.h"

#ifndef CEPETRA_VECTOR_CPP_HPP
#define CEPETRA_VECTOR_CPP_HPP


#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const RCP<Epetra_Vector>
getVector( CT_Epetra_Vector_ID_t id );

/* get Epetra_Vector from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Vector>
getVector( CTrilinos_Object_ID_t id );

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const RCP<const Epetra_Vector>
getConstVector( CT_Epetra_Vector_ID_t id );

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Vector>
getConstVector( CTrilinos_Object_ID_t id );

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
storeVector( Epetra_Vector *pobj );

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
storeConstVector( const Epetra_Vector *pobj );

/* dump contents of Epetra_Vector and const Epetra_Vector tables */
void
purgeVectorTables(  );


} // namespace CEpetra


#endif // CEPETRA_VECTOR_CPP_HPP


