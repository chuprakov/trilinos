
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


#ifdef HAVE_CTRILINOS_AZTECOO


#include "CEpetra_MultiVector_Cpp.hpp"
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#include "CAztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTestMaxIters */
Table<AztecOO_StatusTestMaxIters>& tableOfStatusTestMaxIterss()
{
    static Table<AztecOO_StatusTestMaxIters>
        loc_tableOfStatusTestMaxIterss(CT_AztecOO_StatusTestMaxIters_ID, "CT_AztecOO_StatusTestMaxIters_ID", FALSE);
    return loc_tableOfStatusTestMaxIterss;
}

/* table to hold objects of type const AztecOO_StatusTestMaxIters */
Table<const AztecOO_StatusTestMaxIters>& tableOfConstStatusTestMaxIterss()
{
    static Table<const AztecOO_StatusTestMaxIters>
        loc_tableOfConstStatusTestMaxIterss(CT_AztecOO_StatusTestMaxIters_ID, "CT_AztecOO_StatusTestMaxIters_ID", TRUE);
    return loc_tableOfConstStatusTestMaxIterss;
}


} // namespace


//
// Definitions from CAztecOO_StatusTestMaxIters.h
//


extern "C" {


CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstStatusTestMaxIterss(), id);
    } else {
        newid = CTrilinos::cast(tableOfStatusTestMaxIterss(), id);
    }
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(newid);
}

CTrilinos_Object_ID_t AztecOO_StatusTestMaxIters_Abstract ( 
  CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
}

CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Create ( 
  int MaxIters )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tableOfStatusTestMaxIterss().store(
        new AztecOO_StatusTestMaxIters(MaxIters)));
}

void AztecOO_StatusTestMaxIters_Destroy ( 
  CT_AztecOO_StatusTestMaxIters_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstStatusTestMaxIterss().remove(&aid);
    } else {
        tableOfStatusTestMaxIterss().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(aid);
}

boolean AztecOO_StatusTestMaxIters_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTestMaxIters(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_CheckStatus ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    return (CT_AztecOO_StatusType_E_t)(
         CAztecOO::getStatusTestMaxIters(selfID)->CheckStatus(
        CurrentIter, CEpetra::getMultiVector(
        CurrentResVectorID).getRawPtr(), CurrentResNormEst, ((
        SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_GetStatus ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)(
         CAztecOO::getConstStatusTestMaxIters(
        selfID)->GetStatus() );
}

int AztecOO_StatusTestMaxIters_GetMaxIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(
        selfID)->GetMaxIters();
}

int AztecOO_StatusTestMaxIters_GetNumIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(
        selfID)->GetNumIters();
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestMaxIters_Cpp.hpp
//


/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
    return tableOfStatusTestMaxIterss().get(aid);
}

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CTrilinos_Object_ID_t id )
{
    return tableOfStatusTestMaxIterss().get(id);
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
    if (id.is_const) {
        return tableOfConstStatusTestMaxIterss().get(aid);
    } else {
        return tableOfStatusTestMaxIterss().get(aid);
    }
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstStatusTestMaxIterss().get(id);
    } else {
        return tableOfStatusTestMaxIterss().get(id);
    }
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
            tableOfStatusTestMaxIterss().storeShared(pobj));
}

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
            tableOfConstStatusTestMaxIterss().storeShared(pobj));
}

/* dump contents of AztecOO_StatusTestMaxIters and const AztecOO_StatusTestMaxIters tables */
void
CAztecOO::purgeStatusTestMaxItersTables(  )
{
    tableOfStatusTestMaxIterss().purge();
    tableOfConstStatusTestMaxIterss().purge();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


