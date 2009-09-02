#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_OPERATOR_H
#define CEPETRA_OPERATOR_H


#include "CEpetra_MultiVector.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Operator_Cast() and Epetra_Operator_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Operator table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Operator_ID_t Epetra_Operator_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Operator_Abstract ( 
  CT_Epetra_Operator_ID_t id );

/* Original C++ prototype:
   virtual ~Epetra_Operator();
*/
void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );

/* Original C++ prototype:
   virtual int SetUseTranspose(bool UseTranspose) = 0;
*/
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );

/* Original C++ prototype:
   virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
*/
int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
*/
int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   virtual double NormInf() const = 0;
*/
double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual const char * Label() const = 0;
*/
const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual bool UseTranspose() const = 0;
*/
boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual bool HasNormInf() const = 0;
*/
boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Comm & Comm() const = 0;
*/
CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Map & OperatorDomainMap() const = 0;
*/
CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Map & OperatorRangeMap() const = 0;
*/
CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OPERATOR_H */

