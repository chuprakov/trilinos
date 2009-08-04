#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_TIME_H
#define CEPETRA_TIME_H


#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Time_Cast() and Epetra_Time_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Time table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Time_ID_t Epetra_Time_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Time_Abstract ( 
  CT_Epetra_Time_ID_t id );

/* Original C++ prototype:
   Epetra_Time(const Epetra_Comm & Comm);
*/
CT_Epetra_Time_ID_t Epetra_Time_Create ( 
  CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_Time(const Epetra_Time& Time);
*/
CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( 
  CT_Epetra_Time_ID_t TimeID );

/* Original C++ prototype:
   double WallTime(void) const;
*/
double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID );

/* Original C++ prototype:
   void ResetStartTime(void);
*/
void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID );

/* Original C++ prototype:
   double ElapsedTime(void) const;
*/
double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID );

/* Original C++ prototype:
   virtual ~Epetra_Time(void);
*/
void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID );

/* Original C++ prototype:
   Epetra_Time& operator=(const Epetra_Time& src);
*/
void Epetra_Time_Assign ( 
  CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_TIME_H */

