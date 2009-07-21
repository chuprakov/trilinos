#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_OBJECT_H
#define CEPETRA_OBJECT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Object_Cast() and Epetra_Object_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Object table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Object_ID_t Epetra_Object_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Object_Abstract ( 
  CT_Epetra_Object_ID_t id );

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

/* Original C++ prototype:
   virtual int ReportError(const string Message, int ErrorCode) const;
*/
int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char * const Message, 
  int ErrorCode );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OBJECT_H */

