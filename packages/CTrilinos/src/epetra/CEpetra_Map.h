#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_MAP_H
#define CEPETRA_MAP_H


#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Map_ID_t Epetra_Map_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Map_Abstract ( CT_Epetra_Map_ID_t id );

/* Original C++ prototype:
   Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_Map(const Epetra_Map& map);
*/
CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

/* Original C++ prototype:
   virtual ~Epetra_Map(void);
*/
void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );

/* Original C++ prototype:
   Epetra_Map & operator=(const Epetra_Map & map);
*/
void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MAP_H */

