#include "CTrilinos_config.h"

#ifndef CEPETRA_COMM_CPP_HPP
#define CEPETRA_COMM_CPP_HPP


#include "CEpetra_Comm.h"
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_Comm;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Comm from non-const table using CT_Epetra_Comm_ID */
const RCP<Epetra_Comm>
getComm( CT_Epetra_Comm_ID_t id );

/* get Epetra_Comm from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Comm>
getComm( CTrilinos_Object_ID_t id );

/* get const Epetra_Comm from either the const or non-const table
 * using CT_Epetra_Comm_ID */
const RCP<const Epetra_Comm>
getConstComm( CT_Epetra_Comm_ID_t id );

/* get const Epetra_Comm from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Comm>
getConstComm( CTrilinos_Object_ID_t id );

/* store Epetra_Comm in non-const table */
CT_Epetra_Comm_ID_t
storeComm( Epetra_Comm *pobj );

/* store const Epetra_Comm in const table */
CT_Epetra_Comm_ID_t
storeConstComm( const Epetra_Comm *pobj );

/* dump contents of Epetra_Comm and const Epetra_Comm tables */
void
purgeCommTables(  );


} // namespace CEpetra


#endif // CEPETRA_COMM_CPP_HPP


