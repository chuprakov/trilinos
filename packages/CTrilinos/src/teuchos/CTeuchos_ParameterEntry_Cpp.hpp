#include "CTrilinos_config.h"

#ifndef CTEUCHOS_PARAMETERENTRY_CPP_HPP
#define CTEUCHOS_PARAMETERENTRY_CPP_HPP


#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CTeuchos {


using Teuchos::RCP;


/* get Teuchos::ParameterEntry from non-const table using CT_Teuchos_ParameterEntry_ID */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/* get Teuchos::ParameterEntry from non-const table using CTrilinos_Object_ID_t */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CTrilinos_Object_ID_t id );

/* get const Teuchos::ParameterEntry from either the const or non-const table
 * using CT_Teuchos_ParameterEntry_ID */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/* get const Teuchos::ParameterEntry from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CTrilinos_Object_ID_t id );

/* store Teuchos::ParameterEntry in non-const table */
CT_Teuchos_ParameterEntry_ID_t
storeParameterEntry( Teuchos::ParameterEntry *pobj );

/* store const Teuchos::ParameterEntry in const table */
CT_Teuchos_ParameterEntry_ID_t
storeConstParameterEntry( const Teuchos::ParameterEntry *pobj );

/* dump contents of Teuchos::ParameterEntry and const Teuchos::ParameterEntry tables */
void
purgeParameterEntryTables(  );


} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERENTRY_CPP_HPP


