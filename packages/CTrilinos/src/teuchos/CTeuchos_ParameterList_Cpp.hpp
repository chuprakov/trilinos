#include "CTrilinos_config.h"

#ifndef CTEUCHOS_PARAMETERLIST_CPP_HPP
#define CTEUCHOS_PARAMETERLIST_CPP_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CTeuchos {


using Teuchos::RCP;


/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const RCP<Teuchos::ParameterList>
getParameterList( CT_Teuchos_ParameterList_ID_t id );

/* get Teuchos::ParameterList from non-const table using CTrilinos_Object_ID_t */
const RCP<Teuchos::ParameterList>
getParameterList( CTrilinos_Object_ID_t id );

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CT_Teuchos_ParameterList_ID_t id );

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CTrilinos_Object_ID_t id );

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
storeParameterList( Teuchos::ParameterList *pobj );

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
storeConstParameterList( const Teuchos::ParameterList *pobj );

/* dump contents of Teuchos::ParameterList and const Teuchos::ParameterList tables */
void
purgeParameterListTables(  );


} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERLIST_CPP_HPP


