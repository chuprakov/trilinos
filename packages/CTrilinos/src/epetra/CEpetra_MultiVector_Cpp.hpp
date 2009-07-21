#include "CTrilinos_config.h"

#include "CEpetra_MultiVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_MULTIVECTOR_CPP_HPP
#define CEPETRA_MULTIVECTOR_CPP_HPP


class Epetra_MultiVector;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const RCP<Epetra_MultiVector>
getMultiVector( CT_Epetra_MultiVector_ID_t id );

/* get Epetra_MultiVector from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_MultiVector>
getMultiVector( CTrilinos_Object_ID_t id );

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const RCP<const Epetra_MultiVector>
getConstMultiVector( CT_Epetra_MultiVector_ID_t id );

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_MultiVector>
getConstMultiVector( CTrilinos_Object_ID_t id );

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
storeMultiVector( Epetra_MultiVector *pobj );

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
storeConstMultiVector( const Epetra_MultiVector *pobj );

/* dump contents of Epetra_MultiVector and const Epetra_MultiVector tables */
void
purgeMultiVectorTables(  );


} // namespace CEpetra


#endif // CEPETRA_MULTIVECTOR_CPP_HPP


