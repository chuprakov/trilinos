
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
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Map.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Map.h
//


extern "C" {


CT_Epetra_Map_ID_t Epetra_Map_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Map_Generalize ( 
  CT_Epetra_Map_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
}

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        new Epetra_Map(NumGlobalElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        new Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements, 
        IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_Map> map = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        mapID);
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        new Epetra_Map(*map));
}

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID )
{
    Epetra_Map& self = *( CTrilinos::tableRepos().get<Epetra_Map, 
        CT_Epetra_Map_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_Map> map = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        mapID);
    self = *map;
}


} // extern "C"


//
// Definitions from CEpetra_Map_Cpp.hpp
//


/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Map, CT_Epetra_Map_ID_t>(id);
}

/* get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Map, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Map, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeMap( Epetra_Map *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(pobj, false);
}

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
CEpetra::storeConstMap( const Epetra_Map *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(pobj, false);
}



