#include "CTrilinos_config.h"

#include "CEpetra_CompObject.h"
#include "Epetra_CompObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_COMPOBJECT_CPP_HPP
#define CEPETRA_COMPOBJECT_CPP_HPP


class Epetra_CompObject;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_CompObject>
getCompObject( CT_Epetra_CompObject_ID_t id );

CT_Epetra_CompObject_ID_t
storeCompObject( const Epetra_CompObject *pobj );

void
purgeCompObjectTable(  );


} // namespace CEpetra


#endif // CEPETRA_COMPOBJECT_CPP_HPP


