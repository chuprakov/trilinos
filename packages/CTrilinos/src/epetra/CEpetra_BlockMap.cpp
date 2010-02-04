
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_BlockMap.h"
#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_BlockMap.h
//


extern "C" {


CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_BlockMap_Generalize ( 
  CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id);
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(new Epetra_BlockMap(NumGlobalElements, 
        ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, MyGlobalElements, ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(mapID);
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(new Epetra_BlockMap(*map));
}

void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->RemoteIDList(NumIDs, GIDList, 
        PIDList, LIDList);
}

int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->RemoteIDList(NumIDs, GIDList, 
        PIDList, LIDList, SizeList);
}

int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->LID(GID);
}

int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->GID(LID);
}

int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->FindLocalElementID(PointID, 
        *ElementID, *ElementOffset);
}

boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MyGID(GID_in)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MyLID(LID_in)) ? TRUE : FALSE);
}

int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MinAllGID();
}

int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MaxAllGID();
}

int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MinMyGID();
}

int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MaxMyGID();
}

int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MinLID();
}

int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MaxLID();
}

int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->NumGlobalElements();
}

int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->NumMyElements();
}

int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MyGlobalElements(
        MyGlobalElementList);
}

int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->ElementSize();
}

int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->ElementSize(LID);
}

int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->FirstPointInElement(LID);
}

int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->IndexBase();
}

int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->NumGlobalPoints();
}

int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->NumMyPoints();
}

int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MinMyElementSize();
}

int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MaxMyElementSize();
}

int Epetra_BlockMap_MinElementSize ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MinElementSize();
}

int Epetra_BlockMap_MaxElementSize ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MaxElementSize();
}

boolean Epetra_BlockMap_UniqueGIDs ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->UniqueGIDs()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        selfID)->ConstantElementSize()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(MapID);
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->SameAs(*Map)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(MapID);
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->PointSameAs(*Map)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->LinearMap()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        selfID)->DistributedGlobal()) ? TRUE : FALSE);
}

int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->MyGlobalElements();
}

int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->FirstPointInElementList();
}

int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->ElementSizeList();
}

int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->PointToElementList();
}

int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->ElementSizeList(ElementSizeList);
}

int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->FirstPointInElementList(
        FirstPointInElementList);
}

int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->PointToElementList(
        PointToElementList);
}

CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->Comm() ));
}

boolean Epetra_BlockMap_IsOneToOne ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID)->IsOneToOne()) ? TRUE : FALSE);
}

void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    Epetra_BlockMap& self = *( CTrilinos::tableRepos().get<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_BlockMap> map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(mapID);
    self = *map;
}


} // extern "C"


//
// Definitions from CEpetra_BlockMap_Cpp.hpp
//


/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BlockMap, CT_Epetra_BlockMap_ID_t>(id);
}

/* get Epetra_BlockMap from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BlockMap, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, CT_Epetra_BlockMap_ID_t>(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, CT_Epetra_BlockMap_ID_t>(pobj, false);
}

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeConstBlockMap( const Epetra_BlockMap *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, CT_Epetra_BlockMap_ID_t>(pobj, false);
}



