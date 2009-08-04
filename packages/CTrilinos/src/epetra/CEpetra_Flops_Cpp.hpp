#include "CTrilinos_config.h"

#ifndef CEPETRA_FLOPS_CPP_HPP
#define CEPETRA_FLOPS_CPP_HPP


#include "CEpetra_Flops.h"
#include "Epetra_Flops.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_Flops;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Flops from non-const table using CT_Epetra_Flops_ID */
const RCP<Epetra_Flops>
getFlops( CT_Epetra_Flops_ID_t id );

/* get Epetra_Flops from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Flops>
getFlops( CTrilinos_Object_ID_t id );

/* get const Epetra_Flops from either the const or non-const table
 * using CT_Epetra_Flops_ID */
const RCP<const Epetra_Flops>
getConstFlops( CT_Epetra_Flops_ID_t id );

/* get const Epetra_Flops from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Flops>
getConstFlops( CTrilinos_Object_ID_t id );

/* store Epetra_Flops in non-const table */
CT_Epetra_Flops_ID_t
storeFlops( Epetra_Flops *pobj );

/* store const Epetra_Flops in const table */
CT_Epetra_Flops_ID_t
storeConstFlops( const Epetra_Flops *pobj );

/* dump contents of Epetra_Flops and const Epetra_Flops tables */
void
purgeFlopsTables(  );


} // namespace CEpetra


#endif // CEPETRA_FLOPS_CPP_HPP


