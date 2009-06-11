#include "CEpetra_SerialCommData.h"
#include "Epetra_SerialCommData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_SERIALCOMMDATA_CPP_HPP
#define CEPETRA_SERIALCOMMDATA_CPP_HPP


class Epetra_SerialCommData;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_SerialCommData>
getSerialCommData( CT_Epetra_SerialCommData_ID_t id );

CT_Epetra_SerialCommData_ID_t
storeSerialCommData( const Epetra_SerialCommData *pobj );

void
purgeSerialCommDataTable(  );


} // namespace CEpetra


#endif // CEPETRA_SERIALCOMMDATA_CPP_HPP


