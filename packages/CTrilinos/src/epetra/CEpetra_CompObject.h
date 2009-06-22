
/* True C header file! */


#ifndef CEPETRA_COMPOBJECT_H
#define CEPETRA_COMPOBJECT_H


#include "CEpetra_Flops.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_CompObject_ID_t Epetra_CompObject_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  );

CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID );

void Epetra_CompObject_Destroy ( 
  CT_Epetra_CompObject_ID_t * selfID );

void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID );

void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID );

void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID );

double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID );

void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in );

void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in );

void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in );

void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in );

void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMPOBJECT_H */

