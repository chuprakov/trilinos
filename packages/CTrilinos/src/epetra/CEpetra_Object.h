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

/* Original C++ prototype:
   Epetra_Object(int TracebackModeIn = -1, bool set_label = true);
*/
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );

/* Original C++ prototype:
   Epetra_Object(const char * const Label, int TracebackModeIn = -1);
*/
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );

/* Original C++ prototype:
   Epetra_Object(const Epetra_Object& Object);
*/
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );

/* Original C++ prototype:
   virtual ~Epetra_Object();
*/
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

/* Original C++ prototype:
   virtual void SetLabel(const char * const Label);
*/
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );

/* Original C++ prototype:
   virtual const char * Label() const;
*/
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OBJECT_H */

