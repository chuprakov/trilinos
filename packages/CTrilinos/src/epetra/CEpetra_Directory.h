#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_DIRECTORY_H
#define CEPETRA_DIRECTORY_H


#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Directory_Cast() and Epetra_Directory_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Directory table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Directory_ID_t Epetra_Directory_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Directory_Abstract ( 
  CT_Epetra_Directory_ID_t id );

/* Original C++ prototype:
   virtual ~Epetra_Directory();
*/
void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );

/* Original C++ prototype:
   virtual int GetDirectoryEntries( const Epetra_BlockMap& Map, const int NumEntries, const int * GlobalEntries, int * Procs, int * LocalEntries, int * EntrySizes, bool high_rank_sharing_procs=false) const = 0;
*/
int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs );

/* Original C++ prototype:
   virtual bool GIDsAllUniquelyOwned() const = 0;
*/
boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DIRECTORY_H */

