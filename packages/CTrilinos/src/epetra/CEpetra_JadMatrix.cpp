
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
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_JadMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_JadMatrix */
Table<Epetra_JadMatrix>& tableOfJadMatrixs()
{
    static Table<Epetra_JadMatrix>
        loc_tableOfJadMatrixs(CT_Epetra_JadMatrix_ID, "CT_Epetra_JadMatrix_ID", false);
    return loc_tableOfJadMatrixs;
}

/* table to hold objects of type const Epetra_JadMatrix */
Table<const Epetra_JadMatrix>& tableOfConstJadMatrixs()
{
    static Table<const Epetra_JadMatrix>
        loc_tableOfConstJadMatrixs(CT_Epetra_JadMatrix_ID, "CT_Epetra_JadMatrix_ID", true);
    return loc_tableOfConstJadMatrixs;
}


} // namespace


//
// Definitions from CEpetra_JadMatrix.h
//


extern "C" {


CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstJadMatrixs(), id);
    } else {
        newid = CTrilinos::cast(tableOfJadMatrixs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_JadMatrix_Abstract ( 
  CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
}

CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tableOfJadMatrixs().store(new Epetra_JadMatrix(
        *CEpetra::getConstRowMatrix(MatrixID))));
}

void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstJadMatrixs().remove(&aid);
    } else {
        tableOfJadMatrixs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(aid);
}

int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure )
{
    return CEpetra::getJadMatrix(selfID)->UpdateValues(
        *CEpetra::getConstRowMatrix(MatrixID), ((
        CheckStructure) != FALSE ? true : false));
}

int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values, Indices);
}

int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex )
{
    return CEpetra::getJadMatrix(selfID)->ExtractMyEntryView(
        CurEntry, *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyEntryView(
        CurEntry, *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstJadMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstJadMatrix(selfID)->Multiply(
        ((TransA) != FALSE ? true : false), 
        *CEpetra::getConstMultiVector(XID), *CEpetra::getMultiVector(
        YID));
}

int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstJadMatrix(selfID)->Solve(
        ((Upper) != FALSE ? true : false), ((
        Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), 
        *CEpetra::getConstMultiVector(XID), *CEpetra::getMultiVector(
        YID));
}


} // extern "C"


//
// Definitions from CEpetra_JadMatrix_Cpp.hpp
//


/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
    return tableOfJadMatrixs().get(aid);
}

/* get Epetra_JadMatrix from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfJadMatrixs().get(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
    if (id.is_const) {
        return tableOfConstJadMatrixs().get(aid);
    } else {
        return tableOfJadMatrixs().get(aid);
    }
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstJadMatrixs().get(id);
    } else {
        return tableOfJadMatrixs().get(id);
    }
}

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
            tableOfJadMatrixs().storeShared(pobj));
}

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeConstJadMatrix( const Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
            tableOfConstJadMatrixs().storeShared(pobj));
}

/* dump contents of Epetra_JadMatrix and const Epetra_JadMatrix tables */
void
CEpetra::purgeJadMatrixTables(  )
{
    tableOfJadMatrixs().purge();
    tableOfConstJadMatrixs().purge();
}



