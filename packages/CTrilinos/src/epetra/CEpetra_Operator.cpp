#include "CTrilinos_config.h"

#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Operator.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Operator>& tableOfOperators()
{
    static Table<Epetra_Operator> loc_tableOfOperators(CT_Epetra_Operator_ID, "CT_Epetra_Operator_ID");
    return loc_tableOfOperators;
}


} // namespace


//
// Definitions from CEpetra_Operator.h
//


extern "C" {


CT_Epetra_Operator_ID_t Epetra_Operator_Cast ( 
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfOperators(), id);
}

void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID )
{
    tableOfOperators().remove(selfID);
}

int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose )
{
    return CEpetra::getOperator(selfID)->SetUseTranspose(UseTranspose);
}

int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getOperator(selfID)->Apply(*pX, *pY);
}

int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getOperator(selfID)->ApplyInverse(*pX, *pY);
}

double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getOperator(selfID)->NormInf();
}

const char * Epetra_Operator_Label ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getOperator(selfID)->Label();
}

boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getOperator(selfID)->UseTranspose();
}

boolean Epetra_Operator_HasNormInf ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getOperator(selfID)->HasNormInf();
}

CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeComm(
        &( CEpetra::getOperator(selfID)->Comm() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getOperator(selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getOperator(selfID)->OperatorRangeMap() ));
}


} // extern "C"


//
// Definitions from CEpetra_Operator_Cpp.hpp
//


const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CT_Epetra_Operator_ID_t id )
{
    return tableOfOperators().get(id);
}

CT_Epetra_Operator_ID_t
CEpetra::storeOperator( const Epetra_Operator *pobj )
{
    return tableOfOperators().storeCopy(pobj);
}

void
CEpetra::purgeOperatorTable(  )
{
    tableOfOperators().purge();
}



