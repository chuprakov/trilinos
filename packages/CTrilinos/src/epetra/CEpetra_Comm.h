
/* True C header file! */


#ifndef CEPETRA_COMM_H
#define CEPETRA_COMM_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Comm_ID_t Epetra_Comm_Cast(
  CTrilinos_Object_ID_t id );

void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMM_H */

