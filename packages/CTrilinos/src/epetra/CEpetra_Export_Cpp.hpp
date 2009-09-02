#include "CTrilinos_config.h"

#ifndef CEPETRA_EXPORT_CPP_HPP
#define CEPETRA_EXPORT_CPP_HPP


#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const RCP<Epetra_Export>
getExport( CT_Epetra_Export_ID_t id );

/* get Epetra_Export from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Export>
getExport( CTrilinos_Object_ID_t id );

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const RCP<const Epetra_Export>
getConstExport( CT_Epetra_Export_ID_t id );

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Export>
getConstExport( CTrilinos_Object_ID_t id );

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
storeExport( Epetra_Export *pobj );

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
storeConstExport( const Epetra_Export *pobj );

/* dump contents of Epetra_Export and const Epetra_Export tables */
void
purgeExportTables(  );


} // namespace CEpetra


#endif // CEPETRA_EXPORT_CPP_HPP


