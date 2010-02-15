
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
#include "CTrilinos_TableRepos.hpp"


//
// Definitions from CAztecOO_StatusTestMaxIters.h
//


extern "C" {


CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTestMaxIters_Generalize ( 
  CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
}

CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Create ( 
  int MaxIters )
{
    return CAztecOO::storeNewStatusTestMaxIters(new AztecOO_StatusTestMaxIters(
        MaxIters));
}

void AztecOO_StatusTestMaxIters_Destroy ( 
  CT_AztecOO_StatusTestMaxIters_ID_t * selfID )
{
    CAztecOO::removeStatusTestMaxIters(selfID);
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
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTestMaxIters(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_GetStatus ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTestMaxIters(
        selfID)->GetStatus() );
}

int AztecOO_StatusTestMaxIters_GetMaxIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(selfID)->GetMaxIters();
}

int AztecOO_StatusTestMaxIters_GetNumIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(selfID)->GetNumIters();
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestMaxIters_Cpp.hpp
//


/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestMaxIters>(id);
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestMaxIters>(id);
}

/* store AztecOO_StatusTestMaxIters (owned) in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeNewStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, true));
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* remove AztecOO_StatusTestMaxIters from table using CT_AztecOO_StatusTestMaxIters_ID */
void
CAztecOO::removeStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(aid);
}

/* purge AztecOO_StatusTestMaxIters table */
void
CAztecOO::purgeStatusTestMaxIters(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTestMaxIters>();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


