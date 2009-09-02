#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_LINEARPROBLEM_H
#define CEPETRA_LINEARPROBLEM_H


#include "CEpetra_RowMatrix.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Operator.h"
#include "CEpetra_Vector.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_LinearProblem_Cast() and Epetra_LinearProblem_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_LinearProblem table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_LinearProblem_Abstract ( 
  CT_Epetra_LinearProblem_ID_t id );

/* Original C++ prototype:
   Epetra_LinearProblem(void);
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );

/* Original C++ prototype:
   Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/* Original C++ prototype:
   Epetra_LinearProblem(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/* Original C++ prototype:
   Epetra_LinearProblem(const Epetra_LinearProblem& Problem);
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );

/* Original C++ prototype:
   virtual ~Epetra_LinearProblem(void);
*/
void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID );

/* Original C++ prototype:
   int CheckInput() const;
*/
int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   void AssertSymmetric();
*/
void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   void SetPDL(ProblemDifficultyLevel PDL);
*/
void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL );

/* Original C++ prototype:
   void SetOperator(Epetra_RowMatrix * A);
*/
void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID );

/* Original C++ prototype:
   void SetOperator(Epetra_Operator * A);
*/
void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID );

/* Original C++ prototype:
   void SetLHS(Epetra_MultiVector * X);
*/
void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID );

/* Original C++ prototype:
   void SetRHS(Epetra_MultiVector * B);
*/
void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID );

/* Original C++ prototype:
   int LeftScale(const Epetra_Vector & D);
*/
int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/* Original C++ prototype:
   int RightScale(const Epetra_Vector & D);
*/
int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/* Original C++ prototype:
   Epetra_Operator * GetOperator() const;
*/
CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   Epetra_RowMatrix * GetMatrix() const;
*/
CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   Epetra_MultiVector * GetLHS() const;
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   Epetra_MultiVector * GetRHS() const;
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   ProblemDifficultyLevel GetPDL() const;
*/
CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/* Original C++ prototype:
   bool IsOperatorSymmetric() const;
*/
boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_LINEARPROBLEM_H */

