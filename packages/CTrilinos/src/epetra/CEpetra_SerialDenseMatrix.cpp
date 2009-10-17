
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

#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_SerialDenseMatrix */
Table<Epetra_SerialDenseMatrix>& tableOfSerialDenseMatrixs()
{
    static Table<Epetra_SerialDenseMatrix>
        loc_tableOfSerialDenseMatrixs(CT_Epetra_SerialDenseMatrix_ID, "CT_Epetra_SerialDenseMatrix_ID", FALSE);
    return loc_tableOfSerialDenseMatrixs;
}

/* table to hold objects of type const Epetra_SerialDenseMatrix */
Table<const Epetra_SerialDenseMatrix>& tableOfConstSerialDenseMatrixs()
{
    static Table<const Epetra_SerialDenseMatrix>
        loc_tableOfConstSerialDenseMatrixs(CT_Epetra_SerialDenseMatrix_ID, "CT_Epetra_SerialDenseMatrix_ID", TRUE);
    return loc_tableOfConstSerialDenseMatrixs;
}


} // namespace


//
// Definitions from CEpetra_SerialDenseMatrix.h
//


extern "C" {


CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstSerialDenseMatrixs(), id);
    } else {
        newid = CTrilinos::cast(tableOfSerialDenseMatrixs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_SerialDenseMatrix_Abstract ( 
  CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id);
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_Empty ( 
  boolean set_object_label )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tableOfSerialDenseMatrixs().store(
        new Epetra_SerialDenseMatrix(
        ((set_object_label) != FALSE ? true : false))));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create ( 
  int NumRows, int NumCols, boolean set_object_label )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tableOfSerialDenseMatrixs().store(
        new Epetra_SerialDenseMatrix(NumRows, NumCols, ((
        set_object_label) != FALSE ? true : false))));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, double * A_in, int LDA_in, 
  int NumRows, int NumCols, boolean set_object_label )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tableOfSerialDenseMatrixs().store(
        new Epetra_SerialDenseMatrix((Epetra_DataAccess) CV, A_in, 
        LDA_in, NumRows, NumCols, ((
        set_object_label) != FALSE ? true : false))));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Duplicate ( 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tableOfSerialDenseMatrixs().store(
        new Epetra_SerialDenseMatrix(
        *CEpetra::getConstSerialDenseMatrix(SourceID))));
}

void Epetra_SerialDenseMatrix_Destroy ( 
  CT_Epetra_SerialDenseMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstSerialDenseMatrixs().remove(&aid);
    } else {
        tableOfSerialDenseMatrixs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(aid);
}

int Epetra_SerialDenseMatrix_Shape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Shape(NumRows, 
        NumCols);
}

int Epetra_SerialDenseMatrix_Reshape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Reshape(NumRows, 
        NumCols);
}

int Epetra_SerialDenseMatrix_Multiply_Matrix ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_SerialDenseMatrix_ID_t AID, 
  CT_Epetra_SerialDenseMatrix_ID_t BID, double ScalarThis )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Multiply(TransA, 
        TransB, ScalarAB, *CEpetra::getConstSerialDenseMatrix(AID), 
        *CEpetra::getConstSerialDenseMatrix(BID), ScalarThis);
}

int Epetra_SerialDenseMatrix_Multiply_Vector ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean transA, 
  CT_Epetra_SerialDenseMatrix_ID_t xID, 
  CT_Epetra_SerialDenseMatrix_ID_t yID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Multiply(
        ((transA) != FALSE ? true : false), 
        *CEpetra::getConstSerialDenseMatrix(xID), 
        *CEpetra::getSerialDenseMatrix(yID));
}

int Epetra_SerialDenseMatrix_Scale ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, double ScalarA )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Scale(ScalarA);
}

double Epetra_SerialDenseMatrix_NormOne ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->NormOne();
}

double Epetra_SerialDenseMatrix_NormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->NormInf();
}

int Epetra_SerialDenseMatrix_Random ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Random();
}

int Epetra_SerialDenseMatrix_M ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->M();
}

int Epetra_SerialDenseMatrix_N ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->N();
}

double * Epetra_SerialDenseMatrix_A_Const ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->A();
}

double * Epetra_SerialDenseMatrix_A ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->A();
}

int Epetra_SerialDenseMatrix_LDA ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->LDA();
}

CT_Epetra_DataAccess_E_t Epetra_SerialDenseMatrix_CV ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return (CT_Epetra_DataAccess_E_t)(
         CEpetra::getConstSerialDenseMatrix(selfID)->CV() );
}

double Epetra_SerialDenseMatrix_OneNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->OneNorm();
}

double Epetra_SerialDenseMatrix_InfNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->InfNorm();
}

int Epetra_SerialDenseMatrix_SetUseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean UseTranspose_in )
{
    return CEpetra::getSerialDenseMatrix(selfID)->SetUseTranspose(
        ((UseTranspose_in) != FALSE ? true : false));
}

int Epetra_SerialDenseMatrix_Apply ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Apply(
        *CEpetra::getConstSerialDenseMatrix(XID), 
        *CEpetra::getSerialDenseMatrix(YID));
}

int Epetra_SerialDenseMatrix_ApplyInverse ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->ApplyInverse(
        *CEpetra::getConstSerialDenseMatrix(XID), 
        *CEpetra::getSerialDenseMatrix(YID));
}

const char * Epetra_SerialDenseMatrix_Label ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->Label();
}

boolean Epetra_SerialDenseMatrix_UseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return ((CEpetra::getConstSerialDenseMatrix(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

boolean Epetra_SerialDenseMatrix_HasNormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return ((CEpetra::getConstSerialDenseMatrix(
        selfID)->HasNormInf()) ? TRUE : FALSE);
}

int Epetra_SerialDenseMatrix_RowDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->RowDim();
}

int Epetra_SerialDenseMatrix_ColDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->ColDim();
}

void Epetra_SerialDenseMatrix_Assign ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(selfID) );

    self = *CEpetra::getConstSerialDenseMatrix(SourceID);
}

boolean Epetra_SerialDenseMatrix_IsEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID )
{
    const Epetra_SerialDenseMatrix& self = *( CEpetra::getConstSerialDenseMatrix(selfID) );

    return ((self == *CEpetra::getConstSerialDenseMatrix(
        rhsID)) ? TRUE : FALSE);
}

boolean Epetra_SerialDenseMatrix_NotEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID )
{
    const Epetra_SerialDenseMatrix& self = *( CEpetra::getConstSerialDenseMatrix(selfID) );

    return ((self != *CEpetra::getConstSerialDenseMatrix(
        rhsID)) ? TRUE : FALSE);
}

void Epetra_SerialDenseMatrix_AddTo ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(selfID) );

    self += *CEpetra::getConstSerialDenseMatrix(SourceID);
}

void Epetra_SerialDenseMatrix_setElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex, double * value )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(selfID) );

    self(RowIndex, ColIndex) = *value;
}

double Epetra_SerialDenseMatrix_getElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex )
{
    const Epetra_SerialDenseMatrix& self = *( CEpetra::getConstSerialDenseMatrix(selfID) );

    return self(RowIndex, ColIndex);
}

const double * Epetra_SerialDenseMatrix_getColumn ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int ColIndex )
{
    const Epetra_SerialDenseMatrix& self = *( CEpetra::getConstSerialDenseMatrix(selfID) );

    return self[ColIndex];
}


} // extern "C"


//
// Definitions from CEpetra_SerialDenseMatrix_Cpp.hpp
//


/* get Epetra_SerialDenseMatrix from non-const table using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
CEpetra::getSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id);
    return tableOfSerialDenseMatrixs().get(aid);
}

/* get Epetra_SerialDenseMatrix from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
CEpetra::getSerialDenseMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfSerialDenseMatrixs().get(id);
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
CEpetra::getConstSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id);
    if (id.is_const) {
        return tableOfConstSerialDenseMatrixs().get(aid);
    } else {
        return tableOfSerialDenseMatrixs().get(aid);
    }
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
CEpetra::getConstSerialDenseMatrix( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstSerialDenseMatrixs().get(id);
    } else {
        return tableOfSerialDenseMatrixs().get(id);
    }
}

/* store Epetra_SerialDenseMatrix in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
CEpetra::storeSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
            tableOfSerialDenseMatrixs().storeShared(pobj));
}

/* store const Epetra_SerialDenseMatrix in const table */
CT_Epetra_SerialDenseMatrix_ID_t
CEpetra::storeConstSerialDenseMatrix( const Epetra_SerialDenseMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
            tableOfConstSerialDenseMatrixs().storeShared(pobj));
}

/* dump contents of Epetra_SerialDenseMatrix and const Epetra_SerialDenseMatrix tables */
void
CEpetra::purgeSerialDenseMatrixTables(  )
{
    tableOfSerialDenseMatrixs().purge();
    tableOfConstSerialDenseMatrixs().purge();
}



