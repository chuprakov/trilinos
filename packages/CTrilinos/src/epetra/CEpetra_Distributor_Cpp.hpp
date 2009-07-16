#include "CTrilinos_config.h"

#include "CEpetra_Distributor.h"
#include "Epetra_Distributor.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_DISTRIBUTOR_CPP_HPP
#define CEPETRA_DISTRIBUTOR_CPP_HPP


class Epetra_Distributor;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Distributor>
getDistributor( CT_Epetra_Distributor_ID_t id );

CT_Epetra_Distributor_ID_t
storeDistributor( const Epetra_Distributor *pobj );

void
purgeDistributorTable(  );


} // namespace CEpetra


#endif // CEPETRA_DISTRIBUTOR_CPP_HPP


