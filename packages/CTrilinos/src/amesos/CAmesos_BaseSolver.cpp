
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


#ifdef HAVE_CTRILINOS_AMESOS


#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CAmesos_BaseSolver.h"
#include "Amesos_BaseSolver.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Amesos_BaseSolver */
Table<Amesos_BaseSolver>& tableOfBaseSolvers()
{
    static Table<Amesos_BaseSolver>
        loc_tableOfBaseSolvers(CT_Amesos_BaseSolver_ID, "CT_Amesos_BaseSolver_ID", false);
    return loc_tableOfBaseSolvers;
}

/* table to hold objects of type const Amesos_BaseSolver */
Table<const Amesos_BaseSolver>& tableOfConstBaseSolvers()
{
    static Table<const Amesos_BaseSolver>
        loc_tableOfConstBaseSolvers(CT_Amesos_BaseSolver_ID, "CT_Amesos_BaseSolver_ID", true);
    return loc_tableOfConstBaseSolvers;
}


} // namespace


//
// Definitions from CAmesos_BaseSolver.h
//


extern "C" {


CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstBaseSolvers(), id);
    } else {
        newid = CTrilinos::cast(tableOfBaseSolvers(), id);
    }
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(newid);
}

CTrilinos_Object_ID_t Amesos_BaseSolver_Abstract ( 
  CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id);
}

void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstBaseSolvers().remove(&aid);
    } else {
        tableOfBaseSolvers().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(aid);
}

int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->SymbolicFactorization();
}

int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->NumericFactorization();
}

int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->Solve();
}

int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose )
{
    return CAmesos::getBaseSolver(selfID)->SetUseTranspose(
        ((UseTranspose) != FALSE ? true : false));
}

boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID )
{
    return CAmesos::getBaseSolver(selfID)->SetParameters(
        *CTeuchos::getParameterList(ParameterListID));
}

CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstLinearProblem(
        CAmesos::getConstBaseSolver(selfID)->GetProblem());
}

boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->MatrixShapeOK()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CAmesos::getConstBaseSolver(
        selfID)->Comm() ));
}

int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSymbolicFact();
}

int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumNumericFact();
}

int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSolve();
}

void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintStatus();
}

void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintTiming();
}

void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID )
{
    CAmesos::getBaseSolver(selfID)->setParameterList(
        CTeuchos::getParameterList(paramListID));
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->getNonconstParameterList().getRawPtr());
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->unsetParameterList().getRawPtr());
}

void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID )
{
    CAmesos::getConstBaseSolver(selfID)->GetTiming(
        *CTeuchos::getParameterList(TimingParameterListID));
}


} // extern "C"


//
// Definitions from CAmesos_BaseSolver_Cpp.hpp
//


/* get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id);
    return tableOfBaseSolvers().get(aid);
}

/* get Amesos_BaseSolver from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CTrilinos_Object_ID_t id )
{
    return tableOfBaseSolvers().get(id);
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id);
    if (id.is_const) {
        return tableOfConstBaseSolvers().get(aid);
    } else {
        return tableOfBaseSolvers().get(aid);
    }
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstBaseSolvers().get(id);
    } else {
        return tableOfBaseSolvers().get(id);
    }
}

/* store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeBaseSolver( Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
            tableOfBaseSolvers().storeShared(pobj));
}

/* store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeConstBaseSolver( const Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
            tableOfConstBaseSolvers().storeShared(pobj));
}

/* dump contents of Amesos_BaseSolver and const Amesos_BaseSolver tables */
void
CAmesos::purgeBaseSolverTables(  )
{
    tableOfBaseSolvers().purge();
    tableOfConstBaseSolvers().purge();
}



#endif /* HAVE_CTRILINOS_AMESOS */


