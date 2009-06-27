#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_OBJECT_H
#define CEPETRA_OBJECT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Object_ID_t Epetra_Object_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );

CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );

CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );

void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );

const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OBJECT_H */

