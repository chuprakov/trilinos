#include "CTrilinos_config.h"

#include "CEpetra_Flops.h"
#include "Epetra_Flops.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_FLOPS_CPP_HPP
#define CEPETRA_FLOPS_CPP_HPP


class Epetra_Flops;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Flops>
getFlops( CT_Epetra_Flops_ID_t id );

CT_Epetra_Flops_ID_t
storeFlops( const Epetra_Flops *pobj );

void
purgeFlopsTable(  );


} // namespace CEpetra


#endif // CEPETRA_FLOPS_CPP_HPP


