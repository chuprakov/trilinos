
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


#include "CTrilinos_enums.h"
#include "CAztecOO_StatusTest.h"
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


//
// Definitions from CAztecOO_StatusTest.h
//


extern "C" {


CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTest_Generalize ( 
  CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id);
}

void AztecOO_StatusTest_Destroy ( 
  CT_AztecOO_StatusTest_ID_t * selfID )
{
    CAztecOO::removeStatusTest(selfID);
}

boolean AztecOO_StatusTest_ResidualVectorRequired ( 
  CT_AztecOO_StatusTest_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTest(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTest(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTest(
        selfID)->GetStatus() );
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTest_Cpp.hpp
//


/* get AztecOO_StatusTest from non-const table using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get AztecOO_StatusTest from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTest>(id);
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTest>(id);
}

/* store AztecOO_StatusTest (owned) in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeNewStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, true));
}

/* store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, false));
}

/* store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeConstStatusTest( const AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, false));
}

/* remove AztecOO_StatusTest from table using CT_AztecOO_StatusTest_ID */
void
CAztecOO::removeStatusTest( CT_AztecOO_StatusTest_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(aid);
}

/* purge AztecOO_StatusTest table */
void
CAztecOO::purgeStatusTest(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTest>();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


