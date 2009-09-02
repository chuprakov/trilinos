#include "CTrilinos_config.h"

#ifndef CEPETRA_LINEARPROBLEM_CPP_HPP
#define CEPETRA_LINEARPROBLEM_CPP_HPP


#include "Epetra_LinearProblem.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_LinearProblem from non-const table using CT_Epetra_LinearProblem_ID */
const RCP<Epetra_LinearProblem>
getLinearProblem( CT_Epetra_LinearProblem_ID_t id );

/* get Epetra_LinearProblem from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_LinearProblem>
getLinearProblem( CTrilinos_Object_ID_t id );

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CT_Epetra_LinearProblem_ID */
const RCP<const Epetra_LinearProblem>
getConstLinearProblem( CT_Epetra_LinearProblem_ID_t id );

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_LinearProblem>
getConstLinearProblem( CTrilinos_Object_ID_t id );

/* store Epetra_LinearProblem in non-const table */
CT_Epetra_LinearProblem_ID_t
storeLinearProblem( Epetra_LinearProblem *pobj );

/* store const Epetra_LinearProblem in const table */
CT_Epetra_LinearProblem_ID_t
storeConstLinearProblem( const Epetra_LinearProblem *pobj );

/* dump contents of Epetra_LinearProblem and const Epetra_LinearProblem tables */
void
purgeLinearProblemTables(  );


} // namespace CEpetra


#endif // CEPETRA_LINEARPROBLEM_CPP_HPP


