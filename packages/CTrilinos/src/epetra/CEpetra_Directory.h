
/* True C header file! */


#ifndef CEPETRA_DIRECTORY_H
#define CEPETRA_DIRECTORY_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Directory_ID_t Epetra_Directory_Cast(
  CTrilinos_Object_ID_t id );

void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );

boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DIRECTORY_H */

