#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AMESOS



/* True C header file! */


#ifndef CAMESOS_H
#define CAMESOS_H


#include "CAmesos_BaseSolver.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Amesos_Cast() and Amesos_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Amesos table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/* Original C++ prototype:
   Amesos();
*/
CT_Amesos_ID_t Amesos_Create (  );

/* Original C++ prototype:
   ~Amesos();
*/
void Amesos_Destroy ( CT_Amesos_ID_t * selfID );

/* Original C++ prototype:
   Amesos_BaseSolver *Create(const char *ClassType, const Epetra_LinearProblem& LinearProblem );
*/
CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( 
  CT_Amesos_ID_t selfID, const char * ClassType, 
  CT_Epetra_LinearProblem_ID_t LinearProblemID );

/* Original C++ prototype:
   bool Query(const char * ClassType);
*/
boolean Amesos_Query ( 
  CT_Amesos_ID_t selfID, const char * ClassType );

/* Original C++ prototype:
   static Teuchos::ParameterList GetValidParameters();
*/
CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CAMESOS_H */

#endif /* HAVE_CTRILINOS_AMESOS */


