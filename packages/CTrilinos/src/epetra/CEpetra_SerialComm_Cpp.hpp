#include "CTrilinos_config.h"

#ifndef CEPETRA_SERIALCOMM_CPP_HPP
#define CEPETRA_SERIALCOMM_CPP_HPP


#include "CEpetra_SerialComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_SerialComm;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const RCP<Epetra_SerialComm>
getSerialComm( CT_Epetra_SerialComm_ID_t id );

/* get Epetra_SerialComm from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_SerialComm>
getSerialComm( CTrilinos_Object_ID_t id );

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const RCP<const Epetra_SerialComm>
getConstSerialComm( CT_Epetra_SerialComm_ID_t id );

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_SerialComm>
getConstSerialComm( CTrilinos_Object_ID_t id );

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
storeSerialComm( Epetra_SerialComm *pobj );

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
storeConstSerialComm( const Epetra_SerialComm *pobj );

/* dump contents of Epetra_SerialComm and const Epetra_SerialComm tables */
void
purgeSerialCommTables(  );


} // namespace CEpetra


#endif // CEPETRA_SERIALCOMM_CPP_HPP


