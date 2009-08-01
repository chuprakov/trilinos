#include "CTrilinos_config.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_BlockMap.h"
#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_BlockMap */
Table<Epetra_BlockMap>& tableOfBlockMaps()
{
    static Table<Epetra_BlockMap>
        loc_tableOfBlockMaps(CT_Epetra_BlockMap_ID, "CT_Epetra_BlockMap_ID", false);
    return loc_tableOfBlockMaps;
}

/* table to hold objects of type const Epetra_BlockMap */
Table<const Epetra_BlockMap>& tableOfConstBlockMaps()
{
    static Table<const Epetra_BlockMap>
        loc_tableOfConstBlockMaps(CT_Epetra_BlockMap_ID, "CT_Epetra_BlockMap_ID", true);
    return loc_tableOfConstBlockMaps;
}


} // namespace


//
// Definitions from CEpetra_BlockMap.h
//


extern "C" {


CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstBlockMaps(), id);
    } else {
        newid = CTrilinos::cast(tableOfBlockMaps(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_BlockMap_Abstract ( 
  CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id);
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store(new Epetra_BlockMap(
        NumGlobalElements, ElementSize, IndexBase, *CEpetra::getComm(CommID))));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store(new Epetra_BlockMap(
        NumGlobalElements, NumMyElements, ElementSize, IndexBase, *CEpetra::getComm(CommID))));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store(new Epetra_BlockMap(
        NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, *CEpetra::getComm(CommID))));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store(new Epetra_BlockMap(
        NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, *CEpetra::getComm(CommID))));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store(new Epetra_BlockMap(
        *CEpetra::getBlockMap(mapID))));
}

void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstBlockMaps().remove(&aid);
    } else {
        tableOfBlockMaps().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(aid);
}

int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList )
{
    return CEpetra::getConstBlockMap(selfID)->RemoteIDList(
        NumIDs, GIDList, PIDList, LIDList);
}

int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList )
{
    return CEpetra::getConstBlockMap(selfID)->RemoteIDList(
        NumIDs, GIDList, PIDList, LIDList, SizeList);
}

int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID )
{
    return CEpetra::getConstBlockMap(selfID)->LID(GID);
}

int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->GID(LID);
}

int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset )
{
    return CEpetra::getConstBlockMap(selfID)->FindLocalElementID(
        PointID, *ElementID, *ElementOffset);
}

boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in )
{
    return         CEpetra::getConstBlockMap(selfID)->MyGID(GID_in);
}

boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in )
{
    return         CEpetra::getConstBlockMap(selfID)->MyLID(LID_in);
}

int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinAllGID();
}

int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxAllGID();
}

int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinMyGID();
}

int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxMyGID();
}

int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinLID();
}

int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxLID();
}

int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumGlobalElements();
}

int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumMyElements();
}

int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList )
{
    return CEpetra::getConstBlockMap(selfID)->MyGlobalElements(
        MyGlobalElementList);
}

int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSize();
}

int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSize(LID);
}

int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElement(LID);
}

int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->IndexBase();
}

int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumGlobalPoints();
}

int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumMyPoints();
}

int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinMyElementSize();
}

int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxMyElementSize();
}

int Epetra_BlockMap_MinElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinElementSize();
}

int Epetra_BlockMap_MaxElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxElementSize();
}

boolean Epetra_BlockMap_UniqueGIDs ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return         CEpetra::getConstBlockMap(selfID)->UniqueGIDs();
}

boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return         CEpetra::getConstBlockMap(selfID)->ConstantElementSize();
}

boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    return         CEpetra::getConstBlockMap(selfID)->SameAs(
        *CEpetra::getBlockMap(MapID));
}

boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    return         CEpetra::getConstBlockMap(selfID)->PointSameAs(
        *CEpetra::getBlockMap(MapID));
}

boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID )
{
    return         CEpetra::getConstBlockMap(selfID)->LinearMap();
}

boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return         CEpetra::getConstBlockMap(selfID)->DistributedGlobal();
}

int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MyGlobalElements();
}

int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElementList();
}

int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSizeList();
}

int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->PointToElementList();
}

int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSizeList(
        ElementSizeList);
}

int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElementList(
        FirstPointInElementList);
}

int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList )
{
    return CEpetra::getConstBlockMap(selfID)->PointToElementList(
        PointToElementList);
}

CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::storeConstComm(
        &( CEpetra::getConstBlockMap(selfID)->Comm() ));
}

boolean Epetra_BlockMap_IsOneToOne ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return         CEpetra::getConstBlockMap(selfID)->IsOneToOne();
}

void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    Epetra_BlockMap& self = *( CEpetra::getBlockMap(selfID) );

    self = *CEpetra::getBlockMap(mapID);
}


} // extern "C"


//
// Definitions from CEpetra_BlockMap_Cpp.hpp
//


/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id);
    return tableOfBlockMaps().get(aid);
}

/* get Epetra_BlockMap from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CTrilinos_Object_ID_t id )
{
    return tableOfBlockMaps().get(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id);
    if (id.is_const) {
        return tableOfConstBlockMaps().get(aid);
    } else {
        return tableOfBlockMaps().get(aid);
    }
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstBlockMaps().get(id);
    } else {
        return tableOfBlockMaps().get(id);
    }
}

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
            tableOfBlockMaps().storeCopy(pobj));
}

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeConstBlockMap( const Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
            tableOfConstBlockMaps().storeCopy(pobj));
}

/* dump contents of Epetra_BlockMap and const Epetra_BlockMap tables */
void
CEpetra::purgeBlockMapTables(  )
{
    tableOfBlockMaps().purge();
    tableOfConstBlockMaps().purge();
}



