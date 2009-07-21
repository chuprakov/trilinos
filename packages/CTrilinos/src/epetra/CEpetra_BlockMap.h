#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_BLOCKMAP_H
#define CEPETRA_BLOCKMAP_H


#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_BlockMap_Cast() and Epetra_BlockMap_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_BlockMap table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_BlockMap_Abstract ( 
  CT_Epetra_BlockMap_ID_t id );

/* Original C++ prototype:
   Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, const int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm);
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );

/* Original C++ prototype:
   Epetra_BlockMap(const Epetra_BlockMap& map);
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID );

/* Original C++ prototype:
   virtual ~Epetra_BlockMap(void);
*/
void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID );

/* Original C++ prototype:
   int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList) const;
*/
int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList );

/* Original C++ prototype:
   int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, int * SizeList) const;
*/
int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList );

/* Original C++ prototype:
   int LID(int GID) const;
*/
int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID );

/* Original C++ prototype:
   int GID(int LID) const;
*/
int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID );

/* Original C++ prototype:
   int FindLocalElementID(int PointID, int & ElementID, int & ElementOffset) const;
*/
int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset );

/* Original C++ prototype:
   bool MyGID(int GID_in) const;
*/
boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in );

/* Original C++ prototype:
   bool MyLID(int LID_in) const;
*/
boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in );

/* Original C++ prototype:
   int MinAllGID() const;
*/
int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MaxAllGID() const;
*/
int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MinMyGID() const;
*/
int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MaxMyGID() const;
*/
int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MinLID() const;
*/
int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MaxLID() const;
*/
int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalElements() const;
*/
int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int NumMyElements() const;
*/
int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MyGlobalElements(int * MyGlobalElementList) const;
*/
int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList );

/* Original C++ prototype:
   int ElementSize() const;
*/
int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int ElementSize(int LID) const;
*/
int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );

/* Original C++ prototype:
   int FirstPointInElement(int LID) const;
*/
int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );

/* Original C++ prototype:
   int IndexBase() const;
*/
int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalPoints() const;
*/
int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int NumMyPoints() const;
*/
int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MinMyElementSize() const;
*/
int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MaxMyElementSize() const;
*/
int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MinElementSize() const;
*/
int Epetra_BlockMap_MinElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int MaxElementSize() const;
*/
int Epetra_BlockMap_MaxElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   bool UniqueGIDs() const;
*/
boolean Epetra_BlockMap_UniqueGIDs ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   bool ConstantElementSize() const;
*/
boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   bool SameAs(const Epetra_BlockMap & Map) const;
*/
boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/* Original C++ prototype:
   bool PointSameAs(const Epetra_BlockMap & Map) const;
*/
boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/* Original C++ prototype:
   bool LinearMap() const;
*/
boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   bool DistributedGlobal() const;
*/
boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int * MyGlobalElements() const;
*/
int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int * FirstPointInElementList() const;
*/
int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int * ElementSizeList() const;
*/
int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int * PointToElementList() const;
*/
int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   int ElementSizeList(int * ElementSizeList)const;
*/
int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList );

/* Original C++ prototype:
   int FirstPointInElementList(int * FirstPointInElementList)const;
*/
int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList );

/* Original C++ prototype:
   int PointToElementList(int * PointToElementList) const;
*/
int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList );

/* Original C++ prototype:
   const Epetra_Comm & Comm() const;
*/
CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   bool IsOneToOne() const;
*/
boolean Epetra_BlockMap_IsOneToOne ( 
  CT_Epetra_BlockMap_ID_t selfID );

/* Original C++ prototype:
   Epetra_BlockMap & operator=(const Epetra_BlockMap & map);
*/
void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_BLOCKMAP_H */

