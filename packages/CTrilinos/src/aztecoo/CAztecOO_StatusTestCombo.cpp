
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
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTestCombo */
Table<AztecOO_StatusTestCombo>& tableOfStatusTestCombos()
{
    static Table<AztecOO_StatusTestCombo>
        loc_tableOfStatusTestCombos(CT_AztecOO_StatusTestCombo_ID, "CT_AztecOO_StatusTestCombo_ID", FALSE);
    return loc_tableOfStatusTestCombos;
}

/* table to hold objects of type const AztecOO_StatusTestCombo */
Table<const AztecOO_StatusTestCombo>& tableOfConstStatusTestCombos()
{
    static Table<const AztecOO_StatusTestCombo>
        loc_tableOfConstStatusTestCombos(CT_AztecOO_StatusTestCombo_ID, "CT_AztecOO_StatusTestCombo_ID", TRUE);
    return loc_tableOfConstStatusTestCombos;
}


} // namespace


//
// Definitions from CAztecOO_StatusTestCombo.h
//


extern "C" {


CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstStatusTestCombos(), id);
    } else {
        newid = CTrilinos::cast(tableOfStatusTestCombos(), id);
    }
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(newid);
}

CTrilinos_Object_ID_t AztecOO_StatusTestCombo_Abstract ( 
  CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id);
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( 
  CT_ComboType_E_t t )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tableOfStatusTestCombos().store(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t)));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tableOfStatusTestCombos().store(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, 
        *CAztecOO::getStatusTest(aID))));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID, 
  CT_AztecOO_StatusTest_ID_t bID )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tableOfStatusTestCombos().store(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, 
        *CAztecOO::getStatusTest(aID), *CAztecOO::getStatusTest(
        bID))));
}

void AztecOO_StatusTestCombo_Destroy ( 
  CT_AztecOO_StatusTestCombo_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstStatusTestCombos().remove(&aid);
    } else {
        tableOfStatusTestCombos().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(aid);
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
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id);
    return tableOfStatusTestCombos().get(aid);
}

/* get AztecOO_StatusTestCombo from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<AztecOO_StatusTestCombo>
CAztecOO::getStatusTestCombo( CTrilinos_Object_ID_t id )
{
    return tableOfStatusTestCombos().get(id);
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
CAztecOO::getConstStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id);
    if (id.is_const) {
        return tableOfConstStatusTestCombos().get(aid);
    } else {
        return tableOfStatusTestCombos().get(aid);
    }
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
CAztecOO::getConstStatusTestCombo( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstStatusTestCombos().get(id);
    } else {
        return tableOfStatusTestCombos().get(id);
    }
}

/* store AztecOO_StatusTestCombo in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
            tableOfStatusTestCombos().storeShared(pobj));
}

/* store const AztecOO_StatusTestCombo in const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeConstStatusTestCombo( const AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
            tableOfConstStatusTestCombos().storeShared(pobj));
}

/* dump contents of AztecOO_StatusTestCombo and const AztecOO_StatusTestCombo tables */
void
CAztecOO::purgeStatusTestComboTables(  )
{
    tableOfStatusTestCombos().purge();
    tableOfConstStatusTestCombos().purge();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


