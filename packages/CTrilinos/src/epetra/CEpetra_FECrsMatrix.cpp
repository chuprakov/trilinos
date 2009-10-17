
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
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "CEpetra_FECrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_FECrsMatrix */
Table<Epetra_FECrsMatrix>& tableOfFECrsMatrixs()
{
    static Table<Epetra_FECrsMatrix>
        loc_tableOfFECrsMatrixs(CT_Epetra_FECrsMatrix_ID, "CT_Epetra_FECrsMatrix_ID", FALSE);
    return loc_tableOfFECrsMatrixs;
}

/* table to hold objects of type const Epetra_FECrsMatrix */
Table<const Epetra_FECrsMatrix>& tableOfConstFECrsMatrixs()
{
    static Table<const Epetra_FECrsMatrix>
        loc_tableOfConstFECrsMatrixs(CT_Epetra_FECrsMatrix_ID, "CT_Epetra_FECrsMatrix_ID", TRUE);
    return loc_tableOfConstFECrsMatrixs;
}


} // namespace


//
// Definitions from CEpetra_FECrsMatrix.h
//


extern "C" {


CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstFECrsMatrixs(), id);
    } else {
        newid = CTrilinos::cast(tableOfFECrsMatrixs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_FECrsMatrix_Abstract ( 
  CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id);
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int * NumEntriesPerRow, boolean ignoreNonLocalEntries )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *CEpetra::getConstMap(RowMapID), 
        NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false))));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean ignoreNonLocalEntries )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *CEpetra::getConstMap(RowMapID), 
        NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false))));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int * NumEntriesPerRow, 
  boolean ignoreNonLocalEntries )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *CEpetra::getConstMap(RowMapID), 
        *CEpetra::getConstMap(ColMapID), NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false))));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean ignoreNonLocalEntries )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *CEpetra::getConstMap(RowMapID), 
        *CEpetra::getConstMap(ColMapID), NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false))));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID, 
  boolean ignoreNonLocalEntries )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *CEpetra::getConstCrsGraph(GraphID), 
        ((ignoreNonLocalEntries) != FALSE ? true : false))));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Duplicate ( 
  CT_Epetra_FECrsMatrix_ID_t srcID )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store(new Epetra_FECrsMatrix(
        *CEpetra::getConstFECrsMatrix(srcID))));
}

void Epetra_FECrsMatrix_Destroy ( 
  CT_Epetra_FECrsMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstFECrsMatrixs().remove(&aid);
    } else {
        tableOfFECrsMatrixs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(aid);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_InsertGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        numIndices, indices, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        numRows, rows, numCols, cols, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(indicesID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(rowsID), 
        *CEpetra::getConstIntSerialDenseVector(colsID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(indicesID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(rowsID), 
        *CEpetra::getConstIntSerialDenseVector(colsID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(indicesID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(
        *CEpetra::getConstIntSerialDenseVector(rowsID), 
        *CEpetra::getConstIntSerialDenseVector(colsID), 
        *CEpetra::getConstSerialDenseMatrix(valuesID), format);
}

int Epetra_FECrsMatrix_GlobalAssemble ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean callFillComplete )
{
    return CEpetra::getFECrsMatrix(selfID)->GlobalAssemble(
        ((callFillComplete) != FALSE ? true : false));
}

int Epetra_FECrsMatrix_GlobalAssemble_WithMaps ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_Map_ID_t domain_mapID, CT_Epetra_Map_ID_t range_mapID, 
  boolean callFillComplete )
{
    return CEpetra::getFECrsMatrix(selfID)->GlobalAssemble(
        *CEpetra::getConstMap(domain_mapID), *CEpetra::getConstMap(
        range_mapID), ((callFillComplete) != FALSE ? true : false));
}

void Epetra_FECrsMatrix_setIgnoreNonLocalEntries ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean flag )
{
    CEpetra::getFECrsMatrix(selfID)->setIgnoreNonLocalEntries(
        ((flag) != FALSE ? true : false));
}

void Epetra_FECrsMatrix_Assign ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_FECrsMatrix_ID_t srcID )
{
    Epetra_FECrsMatrix& self = *( CEpetra::getFECrsMatrix(selfID) );

    self = *CEpetra::getConstFECrsMatrix(srcID);
}


} // extern "C"


//
// Definitions from CEpetra_FECrsMatrix_Cpp.hpp
//


/* get Epetra_FECrsMatrix from non-const table using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id);
    return tableOfFECrsMatrixs().get(aid);
}

/* get Epetra_FECrsMatrix from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfFECrsMatrixs().get(id);
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id);
    if (id.is_const) {
        return tableOfConstFECrsMatrixs().get(aid);
    } else {
        return tableOfFECrsMatrixs().get(aid);
    }
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstFECrsMatrixs().get(id);
    } else {
        return tableOfFECrsMatrixs().get(id);
    }
}

/* store Epetra_FECrsMatrix in non-const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
            tableOfFECrsMatrixs().storeShared(pobj));
}

/* store const Epetra_FECrsMatrix in const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeConstFECrsMatrix( const Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
            tableOfConstFECrsMatrixs().storeShared(pobj));
}

/* dump contents of Epetra_FECrsMatrix and const Epetra_FECrsMatrix tables */
void
CEpetra::purgeFECrsMatrixTables(  )
{
    tableOfFECrsMatrixs().purge();
    tableOfConstFECrsMatrixs().purge();
}



