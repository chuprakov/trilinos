#include "CTrilinos_config.h"

#ifndef CEPETRA_LAPACK_CPP_HPP
#define CEPETRA_LAPACK_CPP_HPP


#include "Epetra_LAPACK.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_LAPACK from non-const table using CT_Epetra_LAPACK_ID */
const RCP<Epetra_LAPACK>
getLAPACK( CT_Epetra_LAPACK_ID_t id );

/* get Epetra_LAPACK from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_LAPACK>
getLAPACK( CTrilinos_Object_ID_t id );

/* get const Epetra_LAPACK from either the const or non-const table
 * using CT_Epetra_LAPACK_ID */
const RCP<const Epetra_LAPACK>
getConstLAPACK( CT_Epetra_LAPACK_ID_t id );

/* get const Epetra_LAPACK from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_LAPACK>
getConstLAPACK( CTrilinos_Object_ID_t id );

/* store Epetra_LAPACK in non-const table */
CT_Epetra_LAPACK_ID_t
storeLAPACK( Epetra_LAPACK *pobj );

/* store const Epetra_LAPACK in const table */
CT_Epetra_LAPACK_ID_t
storeConstLAPACK( const Epetra_LAPACK *pobj );

/* dump contents of Epetra_LAPACK and const Epetra_LAPACK tables */
void
purgeLAPACKTables(  );


} // namespace CEpetra


#endif // CEPETRA_LAPACK_CPP_HPP


