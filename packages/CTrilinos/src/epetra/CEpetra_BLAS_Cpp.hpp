#include "CTrilinos_config.h"

#ifndef CEPETRA_BLAS_CPP_HPP
#define CEPETRA_BLAS_CPP_HPP


#include "CEpetra_BLAS.h"
#include "Epetra_BLAS.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_BLAS;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_BLAS from non-const table using CT_Epetra_BLAS_ID */
const RCP<Epetra_BLAS>
getBLAS( CT_Epetra_BLAS_ID_t id );

/* get Epetra_BLAS from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_BLAS>
getBLAS( CTrilinos_Object_ID_t id );

/* get const Epetra_BLAS from either the const or non-const table
 * using CT_Epetra_BLAS_ID */
const RCP<const Epetra_BLAS>
getConstBLAS( CT_Epetra_BLAS_ID_t id );

/* get const Epetra_BLAS from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_BLAS>
getConstBLAS( CTrilinos_Object_ID_t id );

/* store Epetra_BLAS in non-const table */
CT_Epetra_BLAS_ID_t
storeBLAS( Epetra_BLAS *pobj );

/* store const Epetra_BLAS in const table */
CT_Epetra_BLAS_ID_t
storeConstBLAS( const Epetra_BLAS *pobj );

/* dump contents of Epetra_BLAS and const Epetra_BLAS tables */
void
purgeBLASTables(  );


} // namespace CEpetra


#endif // CEPETRA_BLAS_CPP_HPP


