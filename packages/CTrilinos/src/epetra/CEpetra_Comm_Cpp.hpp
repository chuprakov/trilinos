#include "CTrilinos_config.h"

#include "CEpetra_Comm.h"
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_COMM_CPP_HPP
#define CEPETRA_COMM_CPP_HPP


class Epetra_Comm;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Comm>
getComm( CT_Epetra_Comm_ID_t id );

CT_Epetra_Comm_ID_t
storeComm( const Epetra_Comm *pobj );

void
purgeCommTable(  );


} // namespace CEpetra


#endif // CEPETRA_COMM_CPP_HPP


