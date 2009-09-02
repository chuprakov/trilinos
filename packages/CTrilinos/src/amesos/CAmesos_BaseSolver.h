#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AMESOS



/* True C header file! */


#ifndef CAMESOS_BASESOLVER_H
#define CAMESOS_BASESOLVER_H


#include "CTeuchos_ParameterList.h"
#include "CEpetra_LinearProblem.h"
#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Amesos_BaseSolver_Cast() and Amesos_BaseSolver_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Amesos_BaseSolver table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Amesos_BaseSolver_Abstract ( 
  CT_Amesos_BaseSolver_ID_t id );

/* Original C++ prototype:
   virtual ~Amesos_BaseSolver();
*/
void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID );

/* Original C++ prototype:
   virtual int SymbolicFactorization() = 0;
*/
int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int NumericFactorization() = 0;
*/
int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int Solve() = 0;
*/
int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int SetUseTranspose(bool UseTranspose) = 0;
*/
int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose );

/* Original C++ prototype:
   virtual bool UseTranspose() const = 0;
*/
boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0;
*/
int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID );

/* Original C++ prototype:
   virtual const Epetra_LinearProblem* GetProblem() const = 0;
*/
CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual bool MatrixShapeOK() const = 0;
*/
boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Comm & Comm() const = 0;
*/
CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int NumSymbolicFact() const = 0;
*/
int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int NumNumericFact() const = 0;
*/
int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual int NumSolve() const = 0;
*/
int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual void PrintStatus() const = 0;
*/
void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual void PrintTiming() const = 0;
*/
void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
*/
void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID );

/* Original C++ prototype:
   virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
*/
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
*/
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/* Original C++ prototype:
   virtual void GetTiming( Teuchos::ParameterList &TimingParameterList ) const;
*/
void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CAMESOS_BASESOLVER_H */

#endif /* HAVE_CTRILINOS_AMESOS */


