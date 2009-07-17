#include "CTrilinos_config.h"

#include "CEpetra_BLAS.h"
#include "Epetra_BLAS.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_BLAS_CPP_HPP
#define CEPETRA_BLAS_CPP_HPP


class Epetra_BLAS;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_BLAS>
getBLAS( CT_Epetra_BLAS_ID_t id );

const RCP<Epetra_BLAS>
getBLAS( CTrilinos_Object_ID_t id );

CT_Epetra_BLAS_ID_t
storeBLAS( const Epetra_BLAS *pobj );

void
purgeBLASTable(  );


} // namespace CEpetra


#endif // CEPETRA_BLAS_CPP_HPP


