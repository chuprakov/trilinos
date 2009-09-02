#include "CTrilinos_config.h"

#ifndef CEPETRA_OPERATOR_CPP_HPP
#define CEPETRA_OPERATOR_CPP_HPP


#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const RCP<Epetra_Operator>
getOperator( CT_Epetra_Operator_ID_t id );

/* get Epetra_Operator from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Operator>
getOperator( CTrilinos_Object_ID_t id );

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const RCP<const Epetra_Operator>
getConstOperator( CT_Epetra_Operator_ID_t id );

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Operator>
getConstOperator( CTrilinos_Object_ID_t id );

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
storeOperator( Epetra_Operator *pobj );

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
storeConstOperator( const Epetra_Operator *pobj );

/* dump contents of Epetra_Operator and const Epetra_Operator tables */
void
purgeOperatorTables(  );


} // namespace CEpetra


#endif // CEPETRA_OPERATOR_CPP_HPP


