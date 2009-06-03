#include "CEpetra_Map.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_MAP_CPP_HPP
#define CEPETRA_MAP_CPP_HPP


class Epetra_Map;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Map>
getMap( CT_Epetra_Map_ID_t id );

CT_Epetra_Map_ID_t
storeMap( const Epetra_Map *pobj );

void
purgeMapTable(  );


} // namespace CEpetra


#endif // CEPETRA_MAP_CPP_HPP


