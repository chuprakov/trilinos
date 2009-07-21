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


/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const RCP<Epetra_Distributor>
getDistributor( CT_Epetra_Distributor_ID_t id );

/* get Epetra_Distributor from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Distributor>
getDistributor( CTrilinos_Object_ID_t id );

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const RCP<const Epetra_Distributor>
getConstDistributor( CT_Epetra_Distributor_ID_t id );

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Distributor>
getConstDistributor( CTrilinos_Object_ID_t id );

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
storeDistributor( Epetra_Distributor *pobj );

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
storeConstDistributor( const Epetra_Distributor *pobj );

/* dump contents of Epetra_Distributor and const Epetra_Distributor tables */
void
purgeDistributorTables(  );


} // namespace CEpetra


#endif // CEPETRA_DISTRIBUTOR_CPP_HPP


