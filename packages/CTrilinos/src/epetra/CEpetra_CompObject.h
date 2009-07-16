#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_COMPOBJECT_H
#define CEPETRA_COMPOBJECT_H


#include "CEpetra_Flops.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_CompObject_ID_t Epetra_CompObject_Cast ( 
  CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   Epetra_CompObject();
*/
CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  );

/* Original C++ prototype:
   Epetra_CompObject(const Epetra_CompObject& Source);
*/
CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID );

/* Original C++ prototype:
   virtual ~Epetra_CompObject();
*/
void Epetra_CompObject_Destroy ( 
  CT_Epetra_CompObject_ID_t * selfID );

/* Original C++ prototype:
   void SetFlopCounter(const Epetra_Flops & FlopCounter_in);
*/
void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID );

/* Original C++ prototype:
   void SetFlopCounter(const Epetra_CompObject & CompObject);
*/
void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID );

/* Original C++ prototype:
   void UnsetFlopCounter();
*/
void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

/* Original C++ prototype:
   Epetra_Flops * GetFlopCounter() const;
*/
CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

/* Original C++ prototype:
   void ResetFlops() const;
*/
void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID );

/* Original C++ prototype:
   double Flops() const;
*/
double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID );

/* Original C++ prototype:
   void UpdateFlops(int Flops_in) const;
*/
void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in );

/* Original C++ prototype:
   void UpdateFlops(long int Flops_in) const;
*/
void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in );

/* Original C++ prototype:
   void UpdateFlops(double Flops_in) const;
*/
void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in );

/* Original C++ prototype:
   void UpdateFlops(float Flops_in) const;
*/
void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in );

/* Original C++ prototype:
   Epetra_CompObject& operator=(const Epetra_CompObject& src);
*/
void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMPOBJECT_H */

