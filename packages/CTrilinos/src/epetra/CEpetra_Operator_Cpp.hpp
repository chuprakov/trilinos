#include "CTrilinos_config.h"

#include "CEpetra_Operator.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_OPERATOR_CPP_HPP
#define CEPETRA_OPERATOR_CPP_HPP


class Epetra_Operator;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Operator>
getOperator( CT_Epetra_Operator_ID_t id );

CT_Epetra_Operator_ID_t
storeOperator( const Epetra_Operator *pobj );

void
purgeOperatorTable(  );


} // namespace CEpetra


#endif // CEPETRA_OPERATOR_CPP_HPP


