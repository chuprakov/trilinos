
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

#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_CrsMatrix.h
//


extern "C" {


CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CrsMatrix_Generalize ( 
  CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id);
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RowMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix((Epetra_DataAccess) CV, 
        *RowMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RowMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix((Epetra_DataAccess) CV, 
        *RowMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        ColMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix((Epetra_DataAccess) CV, 
        *RowMap, *ColMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        ColMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix((Epetra_DataAccess) CV, 
        *RowMap, *ColMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = 
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(GraphID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix((Epetra_DataAccess) CV, 
        *Graph));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<const Epetra_CrsMatrix> Matrix = 
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(MatrixID);
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(new Epetra_CrsMatrix(*Matrix));
}

void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->PutScalar(ScalarConstant);
}

int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Scale(ScalarConstant);
}

int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InsertGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ReplaceGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->SumIntoGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InsertMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ReplaceMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->SumIntoMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<const Epetra_Vector> Diagonal = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        DiagonalID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ReplaceDiagonalValues(*Diagonal);
}

int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->FillComplete(
        ((OptimizeDataStorage) != FALSE ? true : false));
}

int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage )
{
    const Teuchos::RCP<const Epetra_Map> DomainMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        DomainMapID);
    const Teuchos::RCP<const Epetra_Map> RangeMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RangeMapID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->FillComplete(*DomainMap, *RangeMap, 
        ((OptimizeDataStorage) != FALSE ? true : false));
}

int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->OptimizeStorage();
}

int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MakeDataContiguous();
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> Diagonal = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        DiagonalID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractDiagonalCopy(*Diagonal);
}

int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractMyRowView(MyRow, 
        *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumEntries, *Values);
}

int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ExtractMyRowView(MyRow, 
        *NumEntries, *Values);
}

int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    const Teuchos::RCP<Epetra_Vector> y = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        yID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    const Teuchos::RCP<Epetra_Vector> y = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        yID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Multiply1(((TransA) != 
        FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Multiply1(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    const Teuchos::RCP<Epetra_Vector> y = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        yID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InvRowSums(*x);
}

int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InvRowMaxs(*x);
}

int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->LeftScale(*x);
}

int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InvColSums(*x);
}

int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->InvColMaxs(*x);
}

int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RightScale(*x);
}

boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Filled()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(
        selfID)->StorageOptimized()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(
        selfID)->IndicesAreGlobal()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->IndicesAreLocal()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(
        selfID)->IndicesAreContiguous()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->UpperTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NoDiagonal()) ? TRUE : FALSE);
}

double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NormInf();
}

double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NormOne();
}

double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NormFrobenius();
}

int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumGlobalRows();
}

int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumGlobalCols();
}

int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyNonzeros();
}

int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyRows();
}

int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyCols();
}

int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyDiagonals();
}

int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumGlobalEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumAllocatedGlobalEntries(Row);
}

int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MaxNumEntries();
}

int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->GlobalMaxNumEntries();
}

int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumAllocatedMyEntries(Row);
}

int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->IndexBase();
}

boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->StaticGraph()) ? TRUE : FALSE);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Graph() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RowMap() ));
}

int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(newmapID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ReplaceRowMap(*newmap);
}

boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->HaveColMap()) ? TRUE : FALSE);
}

int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(newmapID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ReplaceColMap(*newmap);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ColMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->DomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Import, CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Comm() ));
}

int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->LRID(GRID_in);
}

int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->GRID(LRID_in);
}

int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->LCID(GCID_in);
}

int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MyGRID(GRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MyLRID(LRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MyGCID(GCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MyLCID(LCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->MyGlobalRow(GID)) ? TRUE : FALSE);
}

const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Label();
}

int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->SetUseTranspose(
        ((UseTranspose_in) != FALSE ? true : false));
}

int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->Apply(*X, *Y);
}

int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ApplyInverse(*X, *Y);
}

boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->HasNormInf()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->UseTranspose()) ? TRUE : FALSE);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->OperatorRangeMap() ));
}

int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Import, CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->RowMatrixImporter());
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->ImportMap() ));
}

int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->TransformToLocal();
}

int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_Map> DomainMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        DomainMapID);
    const Teuchos::RCP<const Epetra_Map> RangeMap = 
        CTrilinos::tableRepos().getConst<Epetra_Map, CT_Epetra_Map_ID_t>(
        RangeMapID);
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID)->TransformToLocal(
        DomainMap.getRawPtr(), RangeMap.getRawPtr());
}

void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID )
{
    Epetra_CrsMatrix& self = *( CTrilinos::tableRepos().get<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_CrsMatrix> src = 
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(srcID);
    self = *src;
}

double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc )
{
    const Epetra_CrsMatrix& self = *( 
        CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, 
        CT_Epetra_CrsMatrix_ID_t>(selfID) );

    return self[Loc];
}


} // extern "C"


//
// Definitions from CEpetra_CrsMatrix_Cpp.hpp
//


/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, CT_Epetra_CrsMatrix_ID_t>(id);
}

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, CT_Epetra_CrsMatrix_ID_t>(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, CT_Epetra_CrsMatrix_ID_t>(pobj, false);
}

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeConstCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_CrsMatrix, CT_Epetra_CrsMatrix_ID_t>(pobj, false);
}



