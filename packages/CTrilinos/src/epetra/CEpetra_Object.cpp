
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
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Object.h
//


extern "C" {


CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label )
{
    return CTrilinos::tableRepos().store<Epetra_Object, 
        CT_Epetra_Object_ID_t>(new Epetra_Object(TracebackModeIn, ((
        set_label) != FALSE ? true : false)));
}

CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn )
{
    return CTrilinos::tableRepos().store<Epetra_Object, 
        CT_Epetra_Object_ID_t>(new Epetra_Object(Label, 
        TracebackModeIn));
}

CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID )
{
    return CTrilinos::tableRepos().store<Epetra_Object, 
        CT_Epetra_Object_ID_t>(new Epetra_Object(
        *CEpetra::getConstObject(ObjectID)));
}

void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
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
    return CTrilinos::tableRepos().get<Epetra_Object, CT_Epetra_Object_ID_t>(id);
}

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CT_Epetra_Object_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Object, CT_Epetra_Object_ID_t>(id);
}

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeObject( Epetra_Object *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Object, CT_Epetra_Object_ID_t>(pobj, false);
}

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
CEpetra::storeConstObject( const Epetra_Object *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Object, CT_Epetra_Object_ID_t>(pobj, false);
}



