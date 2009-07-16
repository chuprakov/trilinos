#include "CTrilinos_config.h"

#include "CEpetra_DistObject.h"
#include "Epetra_DistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_DISTOBJECT_CPP_HPP
#define CEPETRA_DISTOBJECT_CPP_HPP


class Epetra_DistObject;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_DistObject>
getDistObject( CT_Epetra_DistObject_ID_t id );

CT_Epetra_DistObject_ID_t
storeDistObject( const Epetra_DistObject *pobj );

void
purgeDistObjectTable(  );


} // namespace CEpetra


#endif // CEPETRA_DISTOBJECT_CPP_HPP


