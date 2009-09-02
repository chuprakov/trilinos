#include "CTrilinos_config.h"

#ifndef CEPETRA_JADMATRIX_CPP_HPP
#define CEPETRA_JADMATRIX_CPP_HPP


#include "Epetra_JadMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const RCP<Epetra_JadMatrix>
getJadMatrix( CT_Epetra_JadMatrix_ID_t id );

/* get Epetra_JadMatrix from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_JadMatrix>
getJadMatrix( CTrilinos_Object_ID_t id );

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const RCP<const Epetra_JadMatrix>
getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id );

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_JadMatrix>
getConstJadMatrix( CTrilinos_Object_ID_t id );

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
storeJadMatrix( Epetra_JadMatrix *pobj );

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
storeConstJadMatrix( const Epetra_JadMatrix *pobj );

/* dump contents of Epetra_JadMatrix and const Epetra_JadMatrix tables */
void
purgeJadMatrixTables(  );


} // namespace CEpetra


#endif // CEPETRA_JADMATRIX_CPP_HPP


