#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_SRCDISTOBJECT_H
#define CEPETRA_SRCDISTOBJECT_H


#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Cast ( 
  CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   virtual ~Epetra_SrcDistObject();
*/
void Epetra_SrcDistObject_Destroy ( 
  CT_Epetra_SrcDistObject_ID_t * selfID );

/* Original C++ prototype:
   virtual const Epetra_BlockMap & Map() const = 0;
*/
CT_Epetra_BlockMap_ID_t Epetra_SrcDistObject_Map ( 
  CT_Epetra_SrcDistObject_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_SRCDISTOBJECT_H */

