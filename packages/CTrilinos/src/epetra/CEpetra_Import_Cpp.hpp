#include "CTrilinos_config.h"

#include "CEpetra_Import.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_IMPORT_CPP_HPP
#define CEPETRA_IMPORT_CPP_HPP


class Epetra_Import;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Import>
getImport( CT_Epetra_Import_ID_t id );

CT_Epetra_Import_ID_t
storeImport( const Epetra_Import *pobj );

void
purgeImportTable(  );


} // namespace CEpetra


#endif // CEPETRA_IMPORT_CPP_HPP


