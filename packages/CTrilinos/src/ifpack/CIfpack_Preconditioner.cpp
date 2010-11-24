
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

#ifdef HAVE_CTRILINOS_IFPACK

#include "CTrilinos_enums.h"
#include "CIfpack_Preconditioner.h"
#include "CIfpack_Preconditioner_Cpp.hpp"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Ifpack_Preconditioner */
Table<Ifpack_Preconditioner>& tableOfPreconditioners()
{
    static Table<Ifpack_Preconditioner> loc_tableOfPreconditioners(CT_Ifpack_Preconditioner_ID);
    return loc_tableOfPreconditioners;
}


} // namespace


//
// Definitions from CIfpack_Preconditioner.h
//


extern "C" {


CT_Ifpack_Preconditioner_ID_t Ifpack_Preconditioner_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(id);
}

CTrilinos_Universal_ID_t Ifpack_Preconditioner_Generalize ( 
  CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id);
}

int Ifpack_Preconditioner_SetParameters ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    return CIfpack::getPreconditioner(selfID)->SetParameters(*List);
}

int Ifpack_Preconditioner_Initialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getPreconditioner(selfID)->Initialize();
}

boolean Ifpack_Preconditioner_IsInitialized ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return ((CIfpack::getConstPreconditioner(
        selfID)->IsInitialized()) ? TRUE : FALSE);
}

int Ifpack_Preconditioner_Compute ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getPreconditioner(selfID)->Compute();
}

boolean Ifpack_Preconditioner_IsComputed ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return ((CIfpack::getConstPreconditioner(
        selfID)->IsComputed()) ? TRUE : FALSE);
}

double Ifpack_Preconditioner_Condest ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->Condest();
}

int Ifpack_Preconditioner_ApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverse(*X, *Y);
}

CT_Epetra_RowMatrix_ID_t Ifpack_Preconditioner_Matrix ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CEpetra::storeConstRowMatrix(&( CIfpack::getConstPreconditioner(
        selfID)->Matrix() ));
}

int Ifpack_Preconditioner_NumInitialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumInitialize();
}

int Ifpack_Preconditioner_NumCompute ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumCompute();
}

int Ifpack_Preconditioner_NumApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumApplyInverse();
}

double Ifpack_Preconditioner_InitializeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->InitializeTime();
}

double Ifpack_Preconditioner_ComputeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ComputeTime();
}

double Ifpack_Preconditioner_ApplyInverseTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverseTime();
}

double Ifpack_Preconditioner_InitializeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->InitializeFlops();
}

double Ifpack_Preconditioner_ComputeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ComputeFlops();
}

double Ifpack_Preconditioner_ApplyInverseFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverseFlops();
}


} // extern "C"


//
// Definitions from CIfpack_Preconditioner_Cpp.hpp
//


/* get Ifpack_Preconditioner from non-const table using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<Ifpack_Preconditioner>
CIfpack::getPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    if (tableOfPreconditioners().isType(id.table))
        return tableOfPreconditioners().get<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get Ifpack_Preconditioner from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Ifpack_Preconditioner>
CIfpack::getPreconditioner( CTrilinos_Universal_ID_t id )
{
    if (tableOfPreconditioners().isType(id.table))
        return tableOfPreconditioners().get<Ifpack_Preconditioner>(id);
    else
        return CTrilinos::TableRepos::get<Ifpack_Preconditioner>(id);
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<const Ifpack_Preconditioner>
CIfpack::getConstPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    if (tableOfPreconditioners().isType(id.table))
        return tableOfPreconditioners().getConst<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Ifpack_Preconditioner>
CIfpack::getConstPreconditioner( CTrilinos_Universal_ID_t id )
{
    if (tableOfPreconditioners().isType(id.table))
        return tableOfPreconditioners().getConst<Ifpack_Preconditioner>(id);
    else
        return CTrilinos::TableRepos::getConst<Ifpack_Preconditioner>(id);
}

/* store Ifpack_Preconditioner (owned) in non-const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storeNewPreconditioner( Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tableOfPreconditioners().store<Ifpack_Preconditioner>(pobj, true));
}

/* store Ifpack_Preconditioner in non-const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storePreconditioner( Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tableOfPreconditioners().store<Ifpack_Preconditioner>(pobj, false));
}

/* store const Ifpack_Preconditioner in const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storeConstPreconditioner( const Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tableOfPreconditioners().store<Ifpack_Preconditioner>(pobj, false));
}

/* remove Ifpack_Preconditioner from table using CT_Ifpack_Preconditioner_ID */
void
CIfpack::removePreconditioner( CT_Ifpack_Preconditioner_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(*id);
    if (tableOfPreconditioners().isType(aid.table))
        tableOfPreconditioners().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(aid);
}

/* remove Ifpack_Preconditioner from table using CTrilinos_Universal_ID_t */
void
CIfpack::removePreconditioner( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfPreconditioners().isType(aid->table))
        tableOfPreconditioners().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Ifpack_Preconditioner table */
void
CIfpack::purgePreconditioner(  )
{
    tableOfPreconditioners().purge();
}

/* store Ifpack_Preconditioner in non-const table */
CTrilinos_Universal_ID_t
CIfpack::aliasPreconditioner( const Teuchos::RCP< Ifpack_Preconditioner > & robj )
{
    return tableOfPreconditioners().alias(robj);
}

/* store const Ifpack_Preconditioner in const table */
CTrilinos_Universal_ID_t
CIfpack::aliasConstPreconditioner( const Teuchos::RCP< const Ifpack_Preconditioner > & robj )
{
    return tableOfPreconditioners().alias(robj);
}


#endif /* HAVE_CTRILINOS_IFPACK */


