#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_FLOPS_H
#define CEPETRA_FLOPS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Flops_Cast() and Epetra_Flops_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Flops table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Flops_ID_t Epetra_Flops_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Flops_Abstract ( 
  CT_Epetra_Flops_ID_t id );

/* Original C++ prototype:
   Epetra_Flops(void);
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );

/* Original C++ prototype:
   Epetra_Flops(const Epetra_Flops& Flops_in);
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );

/* Original C++ prototype:
   double Flops() const;
*/
double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );

/* Original C++ prototype:
   void ResetFlops();
*/
void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );

/* Original C++ prototype:
   virtual ~Epetra_Flops(void);
*/
void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );

/* Original C++ prototype:
   Epetra_Flops& operator=(const Epetra_Flops& src);
*/
void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_FLOPS_H */

