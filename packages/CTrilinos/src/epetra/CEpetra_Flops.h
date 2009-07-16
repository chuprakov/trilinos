#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_FLOPS_H
#define CEPETRA_FLOPS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Flops_ID_t Epetra_Flops_Cast ( CTrilinos_Object_ID_t id );

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

