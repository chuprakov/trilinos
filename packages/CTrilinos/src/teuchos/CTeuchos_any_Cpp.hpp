#include "CTrilinos_config.h"

#ifndef CTEUCHOS_ANY_CPP_HPP
#define CTEUCHOS_ANY_CPP_HPP


#include "Teuchos_any.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CTeuchos {


using Teuchos::RCP;


/* get Teuchos::any from non-const table using CT_Teuchos_any_ID */
const RCP<Teuchos::any>
getany( CT_Teuchos_any_ID_t id );

/* get Teuchos::any from non-const table using CTrilinos_Object_ID_t */
const RCP<Teuchos::any>
getany( CTrilinos_Object_ID_t id );

/* get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const RCP<const Teuchos::any>
getConstany( CT_Teuchos_any_ID_t id );

/* get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Teuchos::any>
getConstany( CTrilinos_Object_ID_t id );

/* store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
storeany( Teuchos::any *pobj );

/* store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
storeConstany( const Teuchos::any *pobj );

/* dump contents of Teuchos::any and const Teuchos::any tables */
void
purgeanyTables(  );


} // namespace CTeuchos


#endif // CTEUCHOS_ANY_CPP_HPP


