#include "CTrilinos_config.h"

#include "CEpetra_RowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_ROWMATRIX_CPP_HPP
#define CEPETRA_ROWMATRIX_CPP_HPP


class Epetra_RowMatrix;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_RowMatrix>
getRowMatrix( CT_Epetra_RowMatrix_ID_t id );

const RCP<Epetra_RowMatrix>
getRowMatrix( CTrilinos_Object_ID_t id );

CT_Epetra_RowMatrix_ID_t
storeRowMatrix( const Epetra_RowMatrix *pobj );

void
purgeRowMatrixTable(  );


} // namespace CEpetra


#endif // CEPETRA_ROWMATRIX_CPP_HPP


