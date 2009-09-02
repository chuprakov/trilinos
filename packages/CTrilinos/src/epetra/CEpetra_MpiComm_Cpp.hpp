#include "CTrilinos_config.h"


#ifdef HAVE_MPI


#ifndef CEPETRA_MPICOMM_CPP_HPP
#define CEPETRA_MPICOMM_CPP_HPP


#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_MpiComm from non-const table using CT_Epetra_MpiComm_ID */
const RCP<Epetra_MpiComm>
getMpiComm( CT_Epetra_MpiComm_ID_t id );

/* get Epetra_MpiComm from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_MpiComm>
getMpiComm( CTrilinos_Object_ID_t id );

/* get const Epetra_MpiComm from either the const or non-const table
 * using CT_Epetra_MpiComm_ID */
const RCP<const Epetra_MpiComm>
getConstMpiComm( CT_Epetra_MpiComm_ID_t id );

/* get const Epetra_MpiComm from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_MpiComm>
getConstMpiComm( CTrilinos_Object_ID_t id );

/* store Epetra_MpiComm in non-const table */
CT_Epetra_MpiComm_ID_t
storeMpiComm( Epetra_MpiComm *pobj );

/* store const Epetra_MpiComm in const table */
CT_Epetra_MpiComm_ID_t
storeConstMpiComm( const Epetra_MpiComm *pobj );

/* dump contents of Epetra_MpiComm and const Epetra_MpiComm tables */
void
purgeMpiCommTables(  );


} // namespace CEpetra


#endif // CEPETRA_MPICOMM_CPP_HPP


#endif /* HAVE_MPI */


