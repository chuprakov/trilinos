
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
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::any */
Table<Teuchos::any>& tableOfanys()
{
    static Table<Teuchos::any>
        loc_tableOfanys(CT_Teuchos_any_ID, "CT_Teuchos_any_ID", FALSE);
    return loc_tableOfanys;
}

/* table to hold objects of type const Teuchos::any */
Table<const Teuchos::any>& tableOfConstanys()
{
    static Table<const Teuchos::any>
        loc_tableOfConstanys(CT_Teuchos_any_ID, "CT_Teuchos_any_ID", TRUE);
    return loc_tableOfConstanys;
}


} // namespace


//
// Definitions from CTeuchos_any.h
//


extern "C" {


CT_Teuchos_any_ID_t Teuchos_any_Create (  )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store(new Teuchos::any()));
}

CT_Teuchos_any_ID_t Teuchos_any_Create_double ( double value )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store(new Teuchos::any(value)));
}

CT_Teuchos_any_ID_t Teuchos_any_Create_int ( int value )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store(new Teuchos::any(value)));
}

CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store(new Teuchos::any(*CTeuchos::getConstany(
        otherID))));
}

void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Teuchos_any_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstanys().remove(&aid);
    } else {
        tableOfanys().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Teuchos_any_ID_t>(aid);
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
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_any_ID_t>(id);
    return tableOfanys().get(aid);
}

/* get Teuchos::any from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Teuchos::any>
CTeuchos::getany( CTrilinos_Object_ID_t id )
{
    return tableOfanys().get(id);
}

/* get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CT_Teuchos_any_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_any_ID_t>(id);
    if (id.is_const) {
        return tableOfConstanys().get(aid);
    } else {
        return tableOfanys().get(aid);
    }
}

/* get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstanys().get(id);
    } else {
        return tableOfanys().get(id);
    }
}

/* store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
CTeuchos::storeany( Teuchos::any *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
            tableOfanys().storeShared(pobj));
}

/* store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
CTeuchos::storeConstany( const Teuchos::any *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
            tableOfConstanys().storeShared(pobj));
}

/* dump contents of Teuchos::any and const Teuchos::any tables */
void
CTeuchos::purgeanyTables(  )
{
    tableOfanys().purge();
    tableOfConstanys().purge();
}



