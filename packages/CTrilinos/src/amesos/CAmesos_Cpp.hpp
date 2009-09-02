#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AMESOS


#ifndef CAMESOS_CPP_HPP
#define CAMESOS_CPP_HPP


#include "Amesos.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CAmesos {


using Teuchos::RCP;


/* get Amesos from non-const table using CT_Amesos_ID */
const RCP<Amesos>
getAmesos( CT_Amesos_ID_t id );

/* get Amesos from non-const table using CTrilinos_Object_ID_t */
const RCP<Amesos>
getAmesos( CTrilinos_Object_ID_t id );

/* get const Amesos from either the const or non-const table
 * using CT_Amesos_ID */
const RCP<const Amesos>
getConstAmesos( CT_Amesos_ID_t id );

/* get const Amesos from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Amesos>
getConstAmesos( CTrilinos_Object_ID_t id );

/* store Amesos in non-const table */
CT_Amesos_ID_t
storeAmesos( Amesos *pobj );

/* store const Amesos in const table */
CT_Amesos_ID_t
storeConstAmesos( const Amesos *pobj );

/* dump contents of Amesos and const Amesos tables */
void
purgeAmesosTables(  );


} // namespace CAmesos


#endif // CAMESOS_CPP_HPP


#endif /* HAVE_CTRILINOS_AMESOS */


