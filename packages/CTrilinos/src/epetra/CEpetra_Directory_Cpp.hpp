#include "CTrilinos_config.h"

#ifndef CEPETRA_DIRECTORY_CPP_HPP
#define CEPETRA_DIRECTORY_CPP_HPP


#include "CEpetra_Directory.h"
#include "Epetra_Directory.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_Directory;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Directory from non-const table using CT_Epetra_Directory_ID */
const RCP<Epetra_Directory>
getDirectory( CT_Epetra_Directory_ID_t id );

/* get Epetra_Directory from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Directory>
getDirectory( CTrilinos_Object_ID_t id );

/* get const Epetra_Directory from either the const or non-const table
 * using CT_Epetra_Directory_ID */
const RCP<const Epetra_Directory>
getConstDirectory( CT_Epetra_Directory_ID_t id );

/* get const Epetra_Directory from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Directory>
getConstDirectory( CTrilinos_Object_ID_t id );

/* store Epetra_Directory in non-const table */
CT_Epetra_Directory_ID_t
storeDirectory( Epetra_Directory *pobj );

/* store const Epetra_Directory in const table */
CT_Epetra_Directory_ID_t
storeConstDirectory( const Epetra_Directory *pobj );

/* dump contents of Epetra_Directory and const Epetra_Directory tables */
void
purgeDirectoryTables(  );


} // namespace CEpetra


#endif // CEPETRA_DIRECTORY_CPP_HPP


