
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

#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LinearProblem.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"


//
// Definitions from CEpetra_LinearProblem.h
//


extern "C" {


CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_LinearProblem_Generalize ( 
  CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id);
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  )
{
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem());
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_Operator> A = CEpetra::getOperator(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> Problem = 
        CEpetra::getConstLinearProblem(ProblemID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(*Problem));
}

void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID )
{
    CEpetra::removeLinearProblem(selfID);
}

int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::getConstLinearProblem(selfID)->CheckInput();
}

void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    CEpetra::getLinearProblem(selfID)->AssertSymmetric();
}

void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL )
{
    CEpetra::getLinearProblem(selfID)->SetPDL((ProblemDifficultyLevel) PDL);
}

void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    CEpetra::getLinearProblem(selfID)->SetOperator(A.getRawPtr());
}

void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID )
{
    const Teuchos::RCP<Epetra_Operator> A = CEpetra::getOperator(AID);
    CEpetra::getLinearProblem(selfID)->SetOperator(A.getRawPtr());
}

void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID )
{
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    CEpetra::getLinearProblem(selfID)->SetLHS(X.getRawPtr());
}

void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    CEpetra::getLinearProblem(selfID)->SetRHS(B.getRawPtr());
}

int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID )
{
    const Teuchos::RCP<const Epetra_Vector> D = CEpetra::getConstVector(DID);
    return CEpetra::getLinearProblem(selfID)->LeftScale(*D);
}

int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID )
{
    const Teuchos::RCP<const Epetra_Vector> D = CEpetra::getConstVector(DID);
    return CEpetra::getLinearProblem(selfID)->RightScale(*D);
}

CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeOperator(CEpetra::getConstLinearProblem(
        selfID)->GetOperator());
}

CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeRowMatrix(CEpetra::getConstLinearProblem(
        selfID)->GetMatrix());
}

CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeMultiVector(CEpetra::getConstLinearProblem(
        selfID)->GetLHS());
}

CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeMultiVector(CEpetra::getConstLinearProblem(
        selfID)->GetRHS());
}

CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return (CT_ProblemDifficultyLevel_E_t)( CEpetra::getConstLinearProblem(
        selfID)->GetPDL() );
}

boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return ((CEpetra::getConstLinearProblem(
        selfID)->IsOperatorSymmetric()) ? TRUE : FALSE);
}


} // extern "C"


//
// Definitions from CEpetra_LinearProblem_Cpp.hpp
//


/* get Epetra_LinearProblem from non-const table using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get Epetra_LinearProblem from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LinearProblem>(id);
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_LinearProblem>(id);
}

/* store Epetra_LinearProblem (owned) in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeNewLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, true));
}

/* store Epetra_LinearProblem in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, false));
}

/* store const Epetra_LinearProblem in const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeConstLinearProblem( const Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, false));
}

/* remove Epetra_LinearProblem from table using CT_Epetra_LinearProblem_ID */
void
CEpetra::removeLinearProblem( CT_Epetra_LinearProblem_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(aid);
}

/* purge Epetra_LinearProblem table */
void
CEpetra::purgeLinearProblem(  )
{
    CTrilinos::tableRepos().purge<Epetra_LinearProblem>();
}



