
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

#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_Flops.h"
#include "Epetra_Flops.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Flops */
Table<Epetra_Flops>& tableOfFlopss()
{
    static Table<Epetra_Flops>
        loc_tableOfFlopss(CT_Epetra_Flops_ID, "CT_Epetra_Flops_ID", FALSE);
    return loc_tableOfFlopss;
}

/* table to hold objects of type const Epetra_Flops */
Table<const Epetra_Flops>& tableOfConstFlopss()
{
    static Table<const Epetra_Flops>
        loc_tableOfConstFlopss(CT_Epetra_Flops_ID, "CT_Epetra_Flops_ID", TRUE);
    return loc_tableOfConstFlopss;
}


} // namespace


//
// Definitions from CEpetra_Flops.h
//


extern "C" {


CT_Epetra_Flops_ID_t Epetra_Flops_Cast ( CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstFlopss(), id);
    } else {
        newid = CTrilinos::cast(tableOfFlopss(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Flops_Abstract ( 
  CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id);
}

CT_Epetra_Flops_ID_t Epetra_Flops_Create (  )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        tableOfFlopss().store(new Epetra_Flops()));
}

CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        tableOfFlopss().store(new Epetra_Flops(
        *CEpetra::getConstFlops(Flops_inID))));
}

void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstFlopss().remove(&aid);
    } else {
        tableOfFlopss().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(aid);
}

double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID )
{
    return CEpetra::getConstFlops(selfID)->Flops();
}

void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID )
{
    CEpetra::getFlops(selfID)->ResetFlops();
}

void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID )
{
    Epetra_Flops& self = *( CEpetra::getFlops(selfID) );

    self = *CEpetra::getConstFlops(srcID);
}


} // extern "C"


//
// Definitions from CEpetra_Flops_Cpp.hpp
//


/* get Epetra_Flops from non-const table using CT_Epetra_Flops_ID */
const Teuchos::RCP<Epetra_Flops>
CEpetra::getFlops( CT_Epetra_Flops_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id);
    return tableOfFlopss().get(aid);
}

/* get Epetra_Flops from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Flops>
CEpetra::getFlops( CTrilinos_Object_ID_t id )
{
    return tableOfFlopss().get(id);
}

/* get const Epetra_Flops from either the const or non-const table
 * using CT_Epetra_Flops_ID */
const Teuchos::RCP<const Epetra_Flops>
CEpetra::getConstFlops( CT_Epetra_Flops_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id);
    if (id.is_const) {
        return tableOfConstFlopss().get(aid);
    } else {
        return tableOfFlopss().get(aid);
    }
}

/* get const Epetra_Flops from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Flops>
CEpetra::getConstFlops( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstFlopss().get(id);
    } else {
        return tableOfFlopss().get(id);
    }
}

/* store Epetra_Flops in non-const table */
CT_Epetra_Flops_ID_t
CEpetra::storeFlops( Epetra_Flops *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
            tableOfFlopss().storeShared(pobj));
}

/* store const Epetra_Flops in const table */
CT_Epetra_Flops_ID_t
CEpetra::storeConstFlops( const Epetra_Flops *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
            tableOfConstFlopss().storeShared(pobj));
}

/* dump contents of Epetra_Flops and const Epetra_Flops tables */
void
CEpetra::purgeFlopsTables(  )
{
    tableOfFlopss().purge();
    tableOfConstFlopss().purge();
}



