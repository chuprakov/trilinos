
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
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "CAztecOO_StatusTest.h"
#include "AztecOO_StatusTest.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTest */
Table<AztecOO_StatusTest>& tableOfStatusTests()
{
    static Table<AztecOO_StatusTest>
        loc_tableOfStatusTests(CT_AztecOO_StatusTest_ID, "CT_AztecOO_StatusTest_ID", FALSE);
    return loc_tableOfStatusTests;
}

/* table to hold objects of type const AztecOO_StatusTest */
Table<const AztecOO_StatusTest>& tableOfConstStatusTests()
{
    static Table<const AztecOO_StatusTest>
        loc_tableOfConstStatusTests(CT_AztecOO_StatusTest_ID, "CT_AztecOO_StatusTest_ID", TRUE);
    return loc_tableOfConstStatusTests;
}


} // namespace


//
// Definitions from CAztecOO_StatusTest.h
//


extern "C" {


CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstStatusTests(), id);
    } else {
        newid = CTrilinos::cast(tableOfStatusTests(), id);
    }
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(newid);
}

CTrilinos_Object_ID_t AztecOO_StatusTest_Abstract ( 
  CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id);
}

void AztecOO_StatusTest_Destroy ( 
  CT_AztecOO_StatusTest_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstStatusTests().remove(&aid);
    } else {
        tableOfStatusTests().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(aid);
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
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTest(
        selfID)->CheckStatus(CurrentIter, CEpetra::getMultiVector(
        CurrentResVectorID).getRawPtr(), CurrentResNormEst, ((
        SolutionUpdated) != FALSE ? true : false)) );
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
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id);
    return tableOfStatusTests().get(aid);
}

/* get AztecOO_StatusTest from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CTrilinos_Object_ID_t id )
{
    return tableOfStatusTests().get(id);
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id);
    if (id.is_const) {
        return tableOfConstStatusTests().get(aid);
    } else {
        return tableOfStatusTests().get(aid);
    }
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstStatusTests().get(id);
    } else {
        return tableOfStatusTests().get(id);
    }
}

/* store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
            tableOfStatusTests().storeShared(pobj));
}

/* store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeConstStatusTest( const AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
            tableOfConstStatusTests().storeShared(pobj));
}

/* dump contents of AztecOO_StatusTest and const AztecOO_StatusTest tables */
void
CAztecOO::purgeStatusTestTables(  )
{
    tableOfStatusTests().purge();
    tableOfConstStatusTests().purge();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


