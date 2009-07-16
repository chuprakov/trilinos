#include "CTrilinos_config.h"

#include "CEpetra_Export.h"
#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_EXPORT_CPP_HPP
#define CEPETRA_EXPORT_CPP_HPP


class Epetra_Export;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Export>
getExport( CT_Epetra_Export_ID_t id );

CT_Epetra_Export_ID_t
storeExport( const Epetra_Export *pobj );

void
purgeExportTable(  );


} // namespace CEpetra


#endif // CEPETRA_EXPORT_CPP_HPP


