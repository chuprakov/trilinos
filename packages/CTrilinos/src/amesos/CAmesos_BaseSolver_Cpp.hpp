#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AMESOS


#ifndef CAMESOS_BASESOLVER_CPP_HPP
#define CAMESOS_BASESOLVER_CPP_HPP


#include "Amesos_BaseSolver.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CAmesos {


using Teuchos::RCP;


/* get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const RCP<Amesos_BaseSolver>
getBaseSolver( CT_Amesos_BaseSolver_ID_t id );

/* get Amesos_BaseSolver from non-const table using CTrilinos_Object_ID_t */
const RCP<Amesos_BaseSolver>
getBaseSolver( CTrilinos_Object_ID_t id );

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const RCP<const Amesos_BaseSolver>
getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id );

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Amesos_BaseSolver>
getConstBaseSolver( CTrilinos_Object_ID_t id );

/* store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
storeBaseSolver( Amesos_BaseSolver *pobj );

/* store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
storeConstBaseSolver( const Amesos_BaseSolver *pobj );

/* dump contents of Amesos_BaseSolver and const Amesos_BaseSolver tables */
void
purgeBaseSolverTables(  );


} // namespace CAmesos


#endif // CAMESOS_BASESOLVER_CPP_HPP


#endif /* HAVE_CTRILINOS_AMESOS */


