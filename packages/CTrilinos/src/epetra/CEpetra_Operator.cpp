
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

#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Operator.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Operator */
Table<Epetra_Operator>& tableOfOperators()
{
    static Table<Epetra_Operator>
        loc_tableOfOperators(CT_Epetra_Operator_ID, "CT_Epetra_Operator_ID", FALSE);
    return loc_tableOfOperators;
}

/* table to hold objects of type const Epetra_Operator */
Table<const Epetra_Operator>& tableOfConstOperators()
{
    static Table<const Epetra_Operator>
        loc_tableOfConstOperators(CT_Epetra_Operator_ID, "CT_Epetra_Operator_ID", TRUE);
    return loc_tableOfConstOperators;
}


} // namespace


//
// Definitions from CEpetra_Operator.h
//


extern "C" {


CT_Epetra_Operator_ID_t Epetra_Operator_Cast ( 
  CTrilinos_Universal_ID_t id )
{
    CTrilinos_Universal_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstOperators(), id);
    } else {
        newid = CTrilinos::cast(tableOfOperators(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(newid);
}

CTrilinos_Universal_ID_t Epetra_Operator_Abstract ( 
  CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id);
}

void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID )
{
    CTrilinos_Universal_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstOperators().remove(&aid);
    } else {
        tableOfOperators().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(aid);
}

int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose )
{
    return CEpetra::getOperator(selfID)->SetUseTranspose(
        ((UseTranspose) != FALSE ? true : false));
}

int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstOperator(selfID)->Apply(
        *CEpetra::getConstMultiVector(XID), *CEpetra::getMultiVector(
        YID));
}

int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstOperator(selfID)->ApplyInverse(
        *CEpetra::getConstMultiVector(XID), *CEpetra::getMultiVector(
        YID));
}

double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getConstOperator(selfID)->NormInf();
}

const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getConstOperator(selfID)->Label();
}

boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return ((CEpetra::getConstOperator(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return ((CEpetra::getConstOperator(
        selfID)->HasNormInf()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstOperator(
        selfID)->Comm() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstOperator(
        selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstOperator(
        selfID)->OperatorRangeMap() ));
}


} // extern "C"


//
// Definitions from CEpetra_Operator_Cpp.hpp
//


/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CT_Epetra_Operator_ID_t id )
{
    CTrilinos_Universal_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id);
    return tableOfOperators().get(aid);
}

/* get Epetra_Operator from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CTrilinos_Universal_ID_t id )
{
    return tableOfOperators().get(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CT_Epetra_Operator_ID_t id )
{
    CTrilinos_Universal_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id);
    if (id.is_const) {
        return tableOfConstOperators().get(aid);
    } else {
        return tableOfOperators().get(aid);
    }
}

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CTrilinos_Universal_ID_t id )
{
    if (id.is_const) {
        return tableOfConstOperators().get(id);
    } else {
        return tableOfOperators().get(id);
    }
}

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeOperator( Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
            tableOfOperators().storeShared(pobj));
}

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
CEpetra::storeConstOperator( const Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
            tableOfConstOperators().storeShared(pobj));
}

/* dump contents of Epetra_Operator and const Epetra_Operator tables */
void
CEpetra::purgeOperatorTables(  )
{
    tableOfOperators().purge();
    tableOfConstOperators().purge();
}



