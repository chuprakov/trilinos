
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

#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_RowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_RowMatrix.h
//


extern "C" {


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_RowMatrix_Generalize ( 
  CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
}

void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->MaxNumEntries();
}

int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> Diagonal = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        DiagonalID);
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->ExtractDiagonalCopy(*Diagonal);
}

int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->InvRowSums(*x);
}

int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->LeftScale(*x);
}

int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = 
        CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->InvColSums(*x);
}

int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = 
        CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(
        xID);
    return CTrilinos::tableRepos().get<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->RightScale(*x);
}

boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->Filled()) ? TRUE : FALSE);
}

double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NormInf();
}

double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NormOne();
}

int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumGlobalNonzeros();
}

int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumGlobalRows();
}

int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumGlobalCols();
}

int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumGlobalDiagonals();
}

int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumMyNonzeros();
}

int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumMyRows();
}

int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumMyCols();
}

int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->NumMyDiagonals();
}

boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->UpperTriangular()) ? TRUE : FALSE);
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Import, CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(selfID)->RowMatrixImporter());
}


} // extern "C"


//
// Definitions from CEpetra_RowMatrix_Cpp.hpp
//


/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_RowMatrix, CT_Epetra_RowMatrix_ID_t>(id);
}

/* get Epetra_RowMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_RowMatrix, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, CT_Epetra_RowMatrix_ID_t>(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_RowMatrix, CT_Epetra_RowMatrix_ID_t>(pobj, false);
}

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeConstRowMatrix( const Epetra_RowMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_RowMatrix, CT_Epetra_RowMatrix_ID_t>(pobj, false);
}



