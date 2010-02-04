
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


#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CAztecOO_Cpp.hpp"
#include "CAztecOO.h"
#include "AztecOO.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CAztecOO.h
//


extern "C" {


CT_AztecOO_ID_t AztecOO_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_Operator> A = 
        CTrilinos::tableRepos().get<Epetra_Operator, CT_Epetra_Operator_ID_t>(
        AID);
    const Teuchos::RCP<Epetra_MultiVector> X = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(new AztecOO(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = 
        CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(new AztecOO(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_AztecOO_ID_t AztecOO_Create_FromLinearProblem ( 
  CT_Epetra_LinearProblem_ID_t LinearProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> LinearProblem = 
        CTrilinos::tableRepos().getConst<Epetra_LinearProblem, 
        CT_Epetra_LinearProblem_ID_t>(LinearProblemID);
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(new AztecOO(
        *LinearProblem));
}

CT_AztecOO_ID_t AztecOO_Create (  )
{
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(
        new AztecOO());
}

CT_AztecOO_ID_t AztecOO_Duplicate ( CT_AztecOO_ID_t SolverID )
{
    const Teuchos::RCP<const AztecOO> Solver = 
        CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(SolverID);
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(new AztecOO(
        *Solver));
}

void AztecOO_Destroy ( CT_AztecOO_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int AztecOO_SetProblem ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_LinearProblem_ID_t probID, 
  boolean call_SetPrecMatrix )
{
    const Teuchos::RCP<const Epetra_LinearProblem> prob = 
        CTrilinos::tableRepos().getConst<Epetra_LinearProblem, 
        CT_Epetra_LinearProblem_ID_t>(probID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetProblem(*prob, ((call_SetPrecMatrix) != 
        FALSE ? true : false));
}

int AztecOO_SetUserOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t UserOperatorID )
{
    const Teuchos::RCP<Epetra_Operator> UserOperator = 
        CTrilinos::tableRepos().get<Epetra_Operator, CT_Epetra_Operator_ID_t>(
        UserOperatorID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetUserOperator(UserOperator.getRawPtr());
}

int AztecOO_SetUserMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t UserMatrixID, 
  boolean call_SetPrecMatrix )
{
    const Teuchos::RCP<Epetra_RowMatrix> UserMatrix = 
        CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(UserMatrixID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetUserMatrix(UserMatrix.getRawPtr(), ((
        call_SetPrecMatrix) != FALSE ? true : false));
}

int AztecOO_SetLHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t XID )
{
    const Teuchos::RCP<Epetra_MultiVector> X = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetLHS(X.getRawPtr());
}

int AztecOO_SetRHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_MultiVector> B = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetRHS(B.getRawPtr());
}

int AztecOO_SetPrecMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t PrecMatrixID )
{
    const Teuchos::RCP<Epetra_RowMatrix> PrecMatrix = 
        CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(PrecMatrixID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetPrecMatrix(PrecMatrix.getRawPtr());
}

int AztecOO_SetPrecOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t PrecOperatorID )
{
    const Teuchos::RCP<Epetra_Operator> PrecOperator = 
        CTrilinos::tableRepos().get<Epetra_Operator, CT_Epetra_Operator_ID_t>(
        PrecOperatorID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetPrecOperator(PrecOperator.getRawPtr());
}

int AztecOO_ConstructPreconditioner ( 
  CT_AztecOO_ID_t selfID, double * condest )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->ConstructPreconditioner(*condest);
}

int AztecOO_DestroyPreconditioner ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->DestroyPreconditioner();
}

double AztecOO_Condest ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->Condest();
}

int AztecOO_CheckInput ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->CheckInput();
}

CT_Epetra_LinearProblem_ID_t AztecOO_GetProblem ( 
  CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_LinearProblem, 
        CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetProblem());
}

CT_Epetra_Operator_ID_t AztecOO_GetUserOperator ( 
  CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetUserOperator());
}

CT_Epetra_RowMatrix_ID_t AztecOO_GetUserMatrix ( 
  CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetUserMatrix());
}

CT_Epetra_Operator_ID_t AztecOO_GetPrecOperator ( 
  CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetPrecOperator());
}

CT_Epetra_RowMatrix_ID_t AztecOO_GetPrecMatrix ( 
  CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetPrecMatrix());
}

CT_Epetra_MultiVector_ID_t AztecOO_GetLHS ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetLHS());
}

CT_Epetra_MultiVector_ID_t AztecOO_GetRHS ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(CTrilinos::tableRepos().getConst<AztecOO, 
        CT_AztecOO_ID_t>(selfID)->GetRHS());
}

void AztecOO_PrintLinearSystem ( 
  CT_AztecOO_ID_t selfID, const char * name )
{
    CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->PrintLinearSystem(name);
}

int AztecOO_SetParameters ( 
  CT_AztecOO_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t parameterlistID, 
  boolean cerr_warning_if_unused )
{
    const Teuchos::RCP<Teuchos::ParameterList> parameterlist = 
        CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(parameterlistID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetParameters(*parameterlist, ((cerr_warning_if_unused) != 
        FALSE ? true : false));
}

int AztecOO_SetAztecDefaults ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAztecDefaults();
}

int AztecOO_SetAztecOption ( 
  CT_AztecOO_ID_t selfID, int option, int value )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAztecOption(option, value);
}

int AztecOO_GetAztecOption ( CT_AztecOO_ID_t selfID, int option )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetAztecOption(option);
}

int AztecOO_SetAztecParam ( 
  CT_AztecOO_ID_t selfID, int param, double value )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAztecParam(param, value);
}

const int * AztecOO_GetAllAztecOptions ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetAllAztecOptions();
}

const double * AztecOO_GetAllAztecParams ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetAllAztecParams();
}

int AztecOO_SetAllAztecOptions ( 
  CT_AztecOO_ID_t selfID, const int * options )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAllAztecOptions(options);
}

int AztecOO_SetAllAztecParams ( 
  CT_AztecOO_ID_t selfID, const double * params )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAllAztecParams(params);
}

int AztecOO_Iterate_Current ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->Iterate(MaxIters, Tolerance);
}

int AztecOO_Iterate ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID, 
  int MaxIters, double Tolerance )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = 
        CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->Iterate(A.getRawPtr(), X.getRawPtr(), B.getRawPtr(), 
        MaxIters, Tolerance);
}

int AztecOO_recursiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->recursiveIterate(MaxIters, Tolerance);
}

const double * AztecOO_GetAztecStatus ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetAztecStatus();
}

int AztecOO_SetUseAdaptiveDefaultsTrue ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetUseAdaptiveDefaultsTrue();
}

int AztecOO_SetAdaptiveParams ( 
  CT_AztecOO_ID_t selfID, int NumTrials, double * athresholds, 
  double * rthresholds, double condestThreshold, double maxFill, 
  int maxKspace )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SetAdaptiveParams(NumTrials, athresholds, rthresholds, 
        condestThreshold, maxFill, maxKspace);
}

int AztecOO_AdaptiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, int MaxSolveAttempts, 
  double Tolerance )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->AdaptiveIterate(MaxIters, MaxSolveAttempts, Tolerance);
}

int AztecOO_NumIters ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->NumIters();
}

double AztecOO_TrueResidual ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->TrueResidual();
}

double AztecOO_ScaledResidual ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->ScaledResidual();
}

double AztecOO_RecursiveResidual ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->RecursiveResidual();
}

double AztecOO_SolveTime ( CT_AztecOO_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(
        selfID)->SolveTime();
}

int AztecOO_GetAllAztecStatus ( 
  CT_AztecOO_ID_t selfID, double * status )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(
        selfID)->GetAllAztecStatus(status);
}


} // extern "C"


//
// Definitions from CAztecOO_Cpp.hpp
//


/* get AztecOO from non-const table using CT_AztecOO_ID */
const Teuchos::RCP<AztecOO>
CAztecOO::getAztecOO( CT_AztecOO_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO, CT_AztecOO_ID_t>(id);
}

/* get AztecOO from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO>
CAztecOO::getAztecOO( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO, CTrilinos_Universal_ID_t>(id);
}

/* get const AztecOO from either the const or non-const table
 * using CT_AztecOO_ID */
const Teuchos::RCP<const AztecOO>
CAztecOO::getConstAztecOO( CT_AztecOO_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CT_AztecOO_ID_t>(id);
}

/* get const AztecOO from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO>
CAztecOO::getConstAztecOO( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO, CTrilinos_Universal_ID_t>(id);
}

/* store AztecOO in non-const table */
CT_AztecOO_ID_t
CAztecOO::storeAztecOO( AztecOO *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(pobj, false);
}

/* store const AztecOO in const table */
CT_AztecOO_ID_t
CAztecOO::storeConstAztecOO( const AztecOO *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO, CT_AztecOO_ID_t>(pobj, false);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


