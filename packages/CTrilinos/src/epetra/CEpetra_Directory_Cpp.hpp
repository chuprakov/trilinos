#include "CTrilinos_config.h"

#include "CEpetra_Directory.h"
#include "Epetra_Directory.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_DIRECTORY_CPP_HPP
#define CEPETRA_DIRECTORY_CPP_HPP


class Epetra_Directory;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Directory>
getDirectory( CT_Epetra_Directory_ID_t id );

CT_Epetra_Directory_ID_t
storeDirectory( const Epetra_Directory *pobj );

void
purgeDirectoryTable(  );


} // namespace CEpetra


#endif // CEPETRA_DIRECTORY_CPP_HPP


