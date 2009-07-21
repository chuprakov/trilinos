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


/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const RCP<Epetra_RowMatrix>
getRowMatrix( CT_Epetra_RowMatrix_ID_t id );

/* get Epetra_RowMatrix from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_RowMatrix>
getRowMatrix( CTrilinos_Object_ID_t id );

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const RCP<const Epetra_RowMatrix>
getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id );

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_RowMatrix>
getConstRowMatrix( CTrilinos_Object_ID_t id );

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
storeRowMatrix( Epetra_RowMatrix *pobj );

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
storeConstRowMatrix( const Epetra_RowMatrix *pobj );

/* dump contents of Epetra_RowMatrix and const Epetra_RowMatrix tables */
void
purgeRowMatrixTables(  );


} // namespace CEpetra


#endif // CEPETRA_ROWMATRIX_CPP_HPP


