#include "CTrilinos_config.h"

#include "CEpetra_CrsMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_CRSMATRIX_CPP_HPP
#define CEPETRA_CRSMATRIX_CPP_HPP


class Epetra_CrsMatrix;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const RCP<Epetra_CrsMatrix>
getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id );

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_CrsMatrix>
getCrsMatrix( CTrilinos_Object_ID_t id );

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id );

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CTrilinos_Object_ID_t id );

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
storeCrsMatrix( Epetra_CrsMatrix *pobj );

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
storeConstCrsMatrix( const Epetra_CrsMatrix *pobj );

/* dump contents of Epetra_CrsMatrix and const Epetra_CrsMatrix tables */
void
purgeCrsMatrixTables(  );


} // namespace CEpetra


#endif // CEPETRA_CRSMATRIX_CPP_HPP


