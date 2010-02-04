
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
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_JadMatrix.h
//


extern "C" {


CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( 
  CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
}

CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(MatrixID);
    return CTrilinos::tableRepos().store<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(new Epetra_JadMatrix(*Matrix));
}

void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CTrilinos::tableRepos().getConst<Epetra_RowMatrix, 
        CT_Epetra_RowMatrix_ID_t>(MatrixID);
    return CTrilinos::tableRepos().get<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->UpdateValues(*Matrix, ((
        CheckStructure) != FALSE ? true : false));
}

int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex )
{
    return CTrilinos::tableRepos().get<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->ExtractMyEntryView(CurEntry, 
        *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->ExtractMyEntryView(CurEntry, 
        *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, 
        CT_Epetra_JadMatrix_ID_t>(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}


} // extern "C"


//
// Definitions from CEpetra_JadMatrix_Cpp.hpp
//


/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_JadMatrix, CT_Epetra_JadMatrix_ID_t>(id);
}

/* get Epetra_JadMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_JadMatrix, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, CT_Epetra_JadMatrix_ID_t>(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_JadMatrix, CT_Epetra_JadMatrix_ID_t>(pobj, false);
}

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeConstJadMatrix( const Epetra_JadMatrix *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_JadMatrix, CT_Epetra_JadMatrix_ID_t>(pobj, false);
}



