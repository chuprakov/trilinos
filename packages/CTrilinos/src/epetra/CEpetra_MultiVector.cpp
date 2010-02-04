
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

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_MultiVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_MultiVector.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_MultiVector_Generalize ( 
  CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id);
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(MapID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(*Map, NumVectors, ((
        zeroOut) != FALSE ? true : false)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(SourceID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(*Source));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * A, int MyLDA, int NumVectors )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(MapID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Map, A, MyLDA, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(MapID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Map, ArrayOfPointers, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(SourceID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Source, Indices, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(SourceID);
    return CTrilinos::tableRepos().store<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Source, StartIndex, NumVectors));
}

void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReplaceGlobalValue(GlobalRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReplaceGlobalValue(
        GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->SumIntoGlobalValue(GlobalRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->SumIntoGlobalValue(
        GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReplaceMyValue(MyRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReplaceMyValue(MyBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->SumIntoMyValue(MyRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->SumIntoMyValue(MyBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->PutScalar(ScalarConstant);
}

int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Random();
}

int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ExtractCopy(A, MyLDA);
}

int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ExtractCopy(ArrayOfPointers);
}

int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ExtractView(A, MyLDA);
}

int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ExtractView(ArrayOfPointers);
}

int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Dot(*A, Result);
}

int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Abs(*A);
}

int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Reciprocal(*A);
}

int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Scale(ScalarValue);
}

int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Scale(ScalarA, *A);
}

int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Update(ScalarA, *A, ScalarThis);
}

int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Update(ScalarA, *A, ScalarB, *B, 
        ScalarThis);
}

int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Norm1(Result);
}

int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Norm2(Result);
}

int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->NormInf(Result);
}

int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result )
{
    const Teuchos::RCP<const Epetra_MultiVector> Weights = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(WeightsID);
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->NormWeighted(*Weights, Result);
}

int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->MinValue(Result);
}

int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->MaxValue(Result);
}

int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->MeanValue(Result);
}

int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Multiply(TransA, TransB, 
        ScalarAB, *A, *B, ScalarThis);
}

int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Multiply(ScalarAB, *A, *B, 
        ScalarThis);
}

int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(BID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReciprocalMultiply(ScalarAB, *A, 
        *B, ScalarThis);
}

int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->SetSeed(Seed_in);
}

unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Seed();
}

int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->NumVectors();
}

int Epetra_MultiVector_MyLength ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->MyLength();
}

int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->GlobalLength();
}

int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->Stride();
}

boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(
        selfID)->ConstantStride()) ? TRUE : FALSE);
}

int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> map = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(mapID);
    return CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID)->ReplaceMap(*map);
}

void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    Epetra_MultiVector& self = *( 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(SourceID);
    self = *Source;
}

double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID) );

    return self[i];
}

CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(selfID) );

    return CTrilinos::tableRepos().store<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        self(i));
}


} // extern "C"


//
// Definitions from CEpetra_MultiVector_Cpp.hpp
//


/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, CT_Epetra_MultiVector_ID_t>(id);
}

/* get Epetra_MultiVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, CT_Epetra_MultiVector_ID_t>(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_MultiVector, CT_Epetra_MultiVector_ID_t>(pobj, false);
}

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeConstMultiVector( const Epetra_MultiVector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_MultiVector, CT_Epetra_MultiVector_ID_t>(pobj, false);
}



