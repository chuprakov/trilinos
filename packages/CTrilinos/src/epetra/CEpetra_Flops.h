
/* True C header file! */


#ifndef CEPETRA_FLOPS_H
#define CEPETRA_FLOPS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Flops_ID_t Epetra_Flops_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );

CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );

double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );

void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );

void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );

void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_FLOPS_H */

