#include "CEpetra_Object.h"
#include "Epetra_Object.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_OBJECT_CPP_HPP
#define CEPETRA_OBJECT_CPP_HPP


class Epetra_Object;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Object>
getObject( CT_Epetra_Object_ID_t id );

CT_Epetra_Object_ID_t
storeObject( const Epetra_Object *pobj );

void
purgeObjectTable(  );


} // namespace CEpetra


#endif // CEPETRA_OBJECT_CPP_HPP


