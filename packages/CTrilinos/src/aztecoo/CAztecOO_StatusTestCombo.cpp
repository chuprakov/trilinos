
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


#include "CAztecOO_StatusTest_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#include "CAztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestCombo.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CAztecOO_StatusTestCombo.h
//


extern "C" {


CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( 
  CT_ComboType_E_t t )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestCombo, 
        CT_AztecOO_StatusTestCombo_ID_t>(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestCombo, 
        CT_AztecOO_StatusTestCombo_ID_t>(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, 
        *CAztecOO::getStatusTest(aID)));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID, 
  CT_AztecOO_StatusTest_ID_t bID )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestCombo, 
        CT_AztecOO_StatusTestCombo_ID_t>(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, 
        *CAztecOO::getStatusTest(aID), *CAztecOO::getStatusTest(
        bID)));
}

void AztecOO_StatusTestCombo_Destroy ( 
  CT_AztecOO_StatusTestCombo_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_AddStatusTest ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, 
  CT_AztecOO_StatusTest_ID_t aID )
{
    return CAztecOO::storeStatusTestCombo(
        &( CAztecOO::getStatusTestCombo(selfID)->AddStatusTest(
        *CAztecOO::getStatusTest(aID)) ));
}

boolean AztecOO_StatusTestCombo_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTestCombo(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_CheckStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTestCombo(
        selfID)->CheckStatus(CurrentIter, CEpetra::getMultiVector(
        CurrentResVectorID).getRawPtr(), CurrentResNormEst, ((
        SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_GetStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)(
         CAztecOO::getConstStatusTestCombo(selfID)->GetStatus() );
}

CT_ComboType_E_t AztecOO_StatusTestCombo_GetComboType ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return (CT_ComboType_E_t)( CAztecOO::getConstStatusTestCombo(
        selfID)->GetComboType() );
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestCombo_Cpp.hpp
//


/* get AztecOO_StatusTestCombo from non-const table using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<AztecOO_StatusTestCombo>
CAztecOO::getStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestCombo, CT_AztecOO_StatusTestCombo_ID_t>(id);
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
CAztecOO::getConstStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestCombo, CT_AztecOO_StatusTestCombo_ID_t>(id);
}

/* store AztecOO_StatusTestCombo in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestCombo, CT_AztecOO_StatusTestCombo_ID_t>(pobj, false);
}

/* store const AztecOO_StatusTestCombo in const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeConstStatusTestCombo( const AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestCombo, CT_AztecOO_StatusTestCombo_ID_t>(pobj, false);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


