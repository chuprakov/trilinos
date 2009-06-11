#include "CEpetra_CrsMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_CRSMATRIX_CPP_HPP
#define CEPETRA_CRSMATRIX_CPP_HPP


class Epetra_CrsMatrix;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_CrsMatrix>
getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id );

CT_Epetra_CrsMatrix_ID_t
storeCrsMatrix( const Epetra_CrsMatrix *pobj );

void
purgeCrsMatrixTable(  );


} // namespace CEpetra


#endif // CEPETRA_CRSMATRIX_CPP_HPP


