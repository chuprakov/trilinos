
/* True C header file! */


#ifndef CEPETRA_OPERATOR_H
#define CEPETRA_OPERATOR_H


#include "CEpetra_MultiVector.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Operator_ID_t Epetra_Operator_Cast(
  CTrilinos_Object_ID_t id );

void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );

int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );

int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );

const char * Epetra_Operator_Label ( 
  CT_Epetra_Operator_ID_t selfID );

boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );

boolean Epetra_Operator_HasNormInf ( 
  CT_Epetra_Operator_ID_t selfID );

CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );

CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );

CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OPERATOR_H */

