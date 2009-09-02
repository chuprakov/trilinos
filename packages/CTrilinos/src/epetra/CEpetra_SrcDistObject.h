#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_SRCDISTOBJECT_H
#define CEPETRA_SRCDISTOBJECT_H


#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_SrcDistObject_Cast() and Epetra_SrcDistObject_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_SrcDistObject table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_SrcDistObject_Abstract ( 
  CT_Epetra_SrcDistObject_ID_t id );

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

