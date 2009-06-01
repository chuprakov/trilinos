
/* True C header file! */


#ifndef CEPETRA_MAP_H
#define CEPETRA_MAP_H


#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Map_ID_t Epetra_Map_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MAP_H */

