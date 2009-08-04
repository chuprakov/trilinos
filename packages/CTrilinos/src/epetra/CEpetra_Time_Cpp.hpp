#include "CTrilinos_config.h"

#ifndef CEPETRA_TIME_CPP_HPP
#define CEPETRA_TIME_CPP_HPP


#include "CEpetra_Time.h"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_Time;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const RCP<Epetra_Time>
getTime( CT_Epetra_Time_ID_t id );

/* get Epetra_Time from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Time>
getTime( CTrilinos_Object_ID_t id );

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const RCP<const Epetra_Time>
getConstTime( CT_Epetra_Time_ID_t id );

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Time>
getConstTime( CTrilinos_Object_ID_t id );

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
storeTime( Epetra_Time *pobj );

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
storeConstTime( const Epetra_Time *pobj );

/* dump contents of Epetra_Time and const Epetra_Time tables */
void
purgeTimeTables(  );


} // namespace CEpetra


#endif // CEPETRA_TIME_CPP_HPP


