
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
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#include "CAztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestResNorm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CAztecOO_StatusTestResNorm.h
//


extern "C" {


CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTestResNorm_Generalize ( 
  CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Create ( 
  CT_Epetra_Operator_ID_t OperatorID, CT_Epetra_Vector_ID_t LHSID, 
  CT_Epetra_Vector_ID_t RHSID, double Tolerance )
{
    const Teuchos::RCP<const Epetra_Operator> Operator = 
        CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(OperatorID);
    const Teuchos::RCP<const Epetra_Vector> LHS = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        LHSID);
    const Teuchos::RCP<const Epetra_Vector> RHS = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        RHSID);
    return CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(new AztecOO_StatusTestResNorm(
        *Operator, *LHS, *RHS, Tolerance));
}

void AztecOO_StatusTestResNorm_Destroy ( 
  CT_AztecOO_StatusTestResNorm_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int AztecOO_StatusTestResNorm_DefineResForm ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  CT_ResType_E_t TypeOfResidual, CT_NormType_E_t TypeOfNorm, 
  CT_Epetra_Vector_ID_t WeightsID )
{
    const Teuchos::RCP<Epetra_Vector> Weights = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        WeightsID);
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->DefineResForm(
        (AztecOO_StatusTestResNorm::ResType) TypeOfResidual, 
        (AztecOO_StatusTestResNorm::NormType) TypeOfNorm, 
        Weights.getRawPtr());
}

int AztecOO_StatusTestResNorm_DefineScaleForm ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  CT_ScaleType_E_t TypeOfScaling, CT_NormType_E_t TypeOfNorm, 
  CT_Epetra_Vector_ID_t WeightsID, double ScaleValue )
{
    const Teuchos::RCP<Epetra_Vector> Weights = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        WeightsID);
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->DefineScaleForm(
        (AztecOO_StatusTestResNorm::ScaleType) TypeOfScaling, 
        (AztecOO_StatusTestResNorm::NormType) TypeOfNorm, Weights.getRawPtr(), 
        ScaleValue);
}

int AztecOO_StatusTestResNorm_ResetTolerance ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, double Tolerance )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->ResetTolerance(Tolerance);
}

int AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  int maxNumExtraIterations )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->SetMaxNumExtraIterations(
        maxNumExtraIterations);
}

int AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(
        selfID)->GetMaxNumExtraIterations();
}

boolean AztecOO_StatusTestResNorm_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_CheckStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( 
        CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->CheckStatus(CurrentIter, 
        CurrentResVector.getRawPtr(), CurrentResNormEst, ((SolutionUpdated) != 
        FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_GetStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( 
        CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->GetStatus() );
}

void AztecOO_StatusTestResNorm_ResetStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->ResetStatus();
}

double AztecOO_StatusTestResNorm_GetTolerance ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->GetTolerance();
}

double AztecOO_StatusTestResNorm_GetTestValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->GetTestValue();
}

double AztecOO_StatusTestResNorm_GetResNormValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->GetResNormValue();
}

double AztecOO_StatusTestResNorm_GetScaledNormValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, 
        CT_AztecOO_StatusTestResNorm_ID_t>(selfID)->GetScaledNormValue();
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestResNorm_Cpp.hpp
//


/* get AztecOO_StatusTestResNorm from non-const table using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

/* get AztecOO_StatusTestResNorm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm, CTrilinos_Universal_ID_t>(id);
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm, CTrilinos_Universal_ID_t>(id);
}

/* store AztecOO_StatusTestResNorm in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm, CT_AztecOO_StatusTestResNorm_ID_t>(pobj, false);
}

/* store const AztecOO_StatusTestResNorm in const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeConstStatusTestResNorm( const AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm, CT_AztecOO_StatusTestResNorm_ID_t>(pobj, false);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


