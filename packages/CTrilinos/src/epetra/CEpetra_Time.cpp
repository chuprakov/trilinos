
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
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_Time.h"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Time.h
//


extern "C" {


CT_Epetra_Time_ID_t Epetra_Time_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Time_Generalize ( 
  CT_Epetra_Time_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id);
}

CT_Epetra_Time_ID_t Epetra_Time_Create ( CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = 
        CTrilinos::tableRepos().getConst<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        CommID);
    return CTrilinos::tableRepos().store<Epetra_Time, CT_Epetra_Time_ID_t>(
        new Epetra_Time(*Comm));
}

CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( 
  CT_Epetra_Time_ID_t TimeID )
{
    const Teuchos::RCP<const Epetra_Time> Time = 
        CTrilinos::tableRepos().getConst<Epetra_Time, CT_Epetra_Time_ID_t>(
        TimeID);
    return CTrilinos::tableRepos().store<Epetra_Time, CT_Epetra_Time_ID_t>(
        new Epetra_Time(*Time));
}

void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time, CT_Epetra_Time_ID_t>(
        selfID)->WallTime();
}

void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID )
{
    CTrilinos::tableRepos().get<Epetra_Time, CT_Epetra_Time_ID_t>(
        selfID)->ResetStartTime();
}

double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time, CT_Epetra_Time_ID_t>(
        selfID)->ElapsedTime();
}

void Epetra_Time_Assign ( 
  CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID )
{
    Epetra_Time& self = *( CTrilinos::tableRepos().get<Epetra_Time, 
        CT_Epetra_Time_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_Time> src = 
        CTrilinos::tableRepos().getConst<Epetra_Time, CT_Epetra_Time_ID_t>(
        srcID);
    self = *src;
}


} // extern "C"


//
// Definitions from CEpetra_Time_Cpp.hpp
//


/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Time, CT_Epetra_Time_ID_t>(id);
}

/* get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Time, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time, CT_Epetra_Time_ID_t>(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeTime( Epetra_Time *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Time, CT_Epetra_Time_ID_t>(pobj, false);
}

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
CEpetra::storeConstTime( const Epetra_Time *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Time, CT_Epetra_Time_ID_t>(pobj, false);
}



