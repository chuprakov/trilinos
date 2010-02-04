
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

#include "CTeuchos_any_Cpp.hpp"
#include "CTeuchos_any.h"
#include "Teuchos_any.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CTeuchos_any.h
//


extern "C" {


CT_Teuchos_any_ID_t Teuchos_any_Create (  )
{
    return CTrilinos::tableRepos().store<Teuchos::any, 
        CT_Teuchos_any_ID_t>(new Teuchos::any());
}

CT_Teuchos_any_ID_t Teuchos_any_Create_double ( double value )
{
    return CTrilinos::tableRepos().store<Teuchos::any, 
        CT_Teuchos_any_ID_t>(new Teuchos::any(value));
}

CT_Teuchos_any_ID_t Teuchos_any_Create_int ( int value )
{
    return CTrilinos::tableRepos().store<Teuchos::any, 
        CT_Teuchos_any_ID_t>(new Teuchos::any(value));
}

CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID )
{
    return CTrilinos::tableRepos().store<Teuchos::any, 
        CT_Teuchos_any_ID_t>(new Teuchos::any(*CTeuchos::getConstany(
        otherID)));
}

void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID )
{
    return CTeuchos::storeany(&( CTeuchos::getany(selfID)->swap(
        *CTeuchos::getany(rhsID)) ));
}

boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID )
{
    return ((CTeuchos::getConstany(selfID)->empty()) ? TRUE : FALSE);
}

const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID )
{
    return CTeuchos::getConstany(selfID)->typeName().c_str();
}

boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID )
{
    return ((CTeuchos::getConstany(selfID)->same(
        *CTeuchos::getConstany(otherID))) ? TRUE : FALSE);
}

void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID )
{
    Teuchos::any& self = *( CTeuchos::getany(selfID) );

    self = *CTeuchos::getConstany(rhsID);
}


} // extern "C"


//
// Definitions from CTeuchos_any_Cpp.hpp
//


/* get Teuchos::any from non-const table using CT_Teuchos_any_ID */
const Teuchos::RCP<Teuchos::any>
CTeuchos::getany( CT_Teuchos_any_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::any, CT_Teuchos_any_ID_t>(id);
}

/* get Teuchos::any from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::any>
CTeuchos::getany( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::any, CTrilinos_Universal_ID_t>(id);
}

/* get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CT_Teuchos_any_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::any, CT_Teuchos_any_ID_t>(id);
}

/* get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::any, CTrilinos_Universal_ID_t>(id);
}

/* store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
CTeuchos::storeany( Teuchos::any *pobj )
{
    return CTrilinos::tableRepos().store<Teuchos::any, CT_Teuchos_any_ID_t>(pobj, false);
}

/* store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
CTeuchos::storeConstany( const Teuchos::any *pobj )
{
    return CTrilinos::tableRepos().store<Teuchos::any, CT_Teuchos_any_ID_t>(pobj, false);
}



