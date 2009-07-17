#include "CTrilinos_config.h"

#include "CEpetra_SrcDistObject.h"
#include "Epetra_SrcDistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_SRCDISTOBJECT_CPP_HPP
#define CEPETRA_SRCDISTOBJECT_CPP_HPP


class Epetra_SrcDistObject;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_SrcDistObject>
getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id );

const RCP<Epetra_SrcDistObject>
getSrcDistObject( CTrilinos_Object_ID_t id );

CT_Epetra_SrcDistObject_ID_t
storeSrcDistObject( const Epetra_SrcDistObject *pobj );

void
purgeSrcDistObjectTable(  );


} // namespace CEpetra


#endif // CEPETRA_SRCDISTOBJECT_CPP_HPP


