#include "CEpetra_Data.h"
#include "Epetra_Data.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_DATA_CPP_HPP
#define CEPETRA_DATA_CPP_HPP


class Epetra_Data;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Data>
getData( CT_Epetra_Data_ID_t id );

CT_Epetra_Data_ID_t
storeData( const Epetra_Data *pobj );

void
purgeDataTable(  );


} // namespace CEpetra


#endif // CEPETRA_DATA_CPP_HPP


