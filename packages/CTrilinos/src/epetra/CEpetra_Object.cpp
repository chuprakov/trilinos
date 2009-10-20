
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

#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_Object.h"
#include "Epetra_Object.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Object */
Table<Epetra_Object>& tableOfObjects()
{
    static Table<Epetra_Object>
        loc_tableOfObjects(CT_Epetra_Object_ID, "CT_Epetra_Object_ID", FALSE);
    return loc_tableOfObjects;
}

/* table to hold objects of type const Epetra_Object */
Table<const Epetra_Object>& tableOfConstObjects()
{
    static Table<const Epetra_Object>
        loc_tableOfConstObjects(CT_Epetra_Object_ID, "CT_Epetra_Object_ID", TRUE);
    return loc_tableOfConstObjects;
}


} // namespace


//
// Definitions from CEpetra_Object.h
//


extern "C" {


CT_Epetra_Object_ID_t Epetra_Object_Cast ( CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstObjects(), id);
    } else {
        newid = CTrilinos::cast(tableOfObjects(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Object_Abstract ( 
  CT_Epetra_Object_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id);
}

CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store(new Epetra_Object(TracebackModeIn, ((
        set_label) != FALSE ? true : false))));
}

CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store(new Epetra_Object(Label, 
        TracebackModeIn)));
}

CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store(new Epetra_Object(
        *CEpetra::getConstObject(ObjectID))));
}

void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Object_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstObjects().remove(&aid);
    } else {
        tableOfObjects().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Object_ID_t>(aid);
}

void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label )
{
    CEpetra::getObject(selfID)->SetLabel(Label);
}

const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID )
{
    return CEpetra::getConstObject(selfID)->Label();
}

int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char Message[], int ErrorCode )
{
    return CEpetra::getConstObject(selfID)->ReportError(std::string(
        Message), ErrorCode);
}

void Epetra_Object_SetTracebackMode ( int TracebackModeValue )
{
    Epetra_Object::SetTracebackMode(TracebackModeValue);
}

int Epetra_Object_GetTracebackMode (  )
{
    return Epetra_Object::GetTracebackMode();
}


} // extern "C"


//
// Definitions from CEpetra_Object_Cpp.hpp
//


/* get Epetra_Object from non-const table using CT_Epetra_Object_ID */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CT_Epetra_Object_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id);
    return tableOfObjects().get(aid);
}

/* get Epetra_Object from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CTrilinos_Object_ID_t id )
{
    return tableOfObjects().get(id);
}

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CT_Epetra_Object_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id);
    if (id.is_const) {
        return tableOfConstObjects().get(aid);
    } else {
        return tableOfObjects().get(aid);
    }
}

/* get const Epetra_Object from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstObjects().get(id);
    } else {
        return tableOfObjects().get(id);
    }
}

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeObject( Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
            tableOfObjects().storeShared(pobj));
}

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
CEpetra::storeConstObject( const Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
            tableOfConstObjects().storeShared(pobj));
}

/* dump contents of Epetra_Object and const Epetra_Object tables */
void
CEpetra::purgeObjectTables(  )
{
    tableOfObjects().purge();
    tableOfConstObjects().purge();
}



