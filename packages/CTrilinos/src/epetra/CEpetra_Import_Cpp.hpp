#include "CTrilinos_config.h"

#ifndef CEPETRA_IMPORT_CPP_HPP
#define CEPETRA_IMPORT_CPP_HPP


#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Import from non-const table using CT_Epetra_Import_ID */
const RCP<Epetra_Import>
getImport( CT_Epetra_Import_ID_t id );

/* get Epetra_Import from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Import>
getImport( CTrilinos_Object_ID_t id );

/* get const Epetra_Import from either the const or non-const table
 * using CT_Epetra_Import_ID */
const RCP<const Epetra_Import>
getConstImport( CT_Epetra_Import_ID_t id );

/* get const Epetra_Import from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Import>
getConstImport( CTrilinos_Object_ID_t id );

/* store Epetra_Import in non-const table */
CT_Epetra_Import_ID_t
storeImport( Epetra_Import *pobj );

/* store const Epetra_Import in const table */
CT_Epetra_Import_ID_t
storeConstImport( const Epetra_Import *pobj );

/* dump contents of Epetra_Import and const Epetra_Import tables */
void
purgeImportTables(  );


} // namespace CEpetra


#endif // CEPETRA_IMPORT_CPP_HPP


