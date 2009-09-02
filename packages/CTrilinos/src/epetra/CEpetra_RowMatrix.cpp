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
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_RowMatrix */
Table<Epetra_RowMatrix>& tableOfRowMatrixs()
{
    static Table<Epetra_RowMatrix>
        loc_tableOfRowMatrixs(CT_Epetra_RowMatrix_ID, "CT_Epetra_RowMatrix_ID", false);
    return loc_tableOfRowMatrixs;
}

/* table to hold objects of type const Epetra_RowMatrix */
Table<const Epetra_RowMatrix>& tableOfConstRowMatrixs()
{
    static Table<const Epetra_RowMatrix>
        loc_tableOfConstRowMatrixs(CT_Epetra_RowMatrix_ID, "CT_Epetra_RowMatrix_ID", true);
    return loc_tableOfConstRowMatrixs;
}


} // namespace


//
// Definitions from CEpetra_RowMatrix.h
//


extern "C" {


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstRowMatrixs(), id);
    } else {
        newid = CTrilinos::cast(tableOfRowMatrixs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_RowMatrix_Abstract ( 
  CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
}

void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstRowMatrixs().remove(&aid);
    } else {
        tableOfRowMatrixs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(aid);
}

int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->MaxNumEntries();
}

int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstRowMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values, Indices);
}

int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    return CEpetra::getConstRowMatrix(selfID)->ExtractDiagonalCopy(
        *CEpetra::getVector(DiagonalID));
}

int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstRowMatrix(selfID)->Multiply(TransA, 
        *CEpetra::getConstMultiVector(XID), *CEpetra::getMultiVector(
        YID));
}

int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstRowMatrix(selfID)->Solve(Upper, Trans, 
        UnitDiagonal, *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstRowMatrix(selfID)->InvRowSums(
        *CEpetra::getVector(xID));
}

int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getRowMatrix(selfID)->LeftScale(
        *CEpetra::getConstVector(xID));
}

int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstRowMatrix(selfID)->InvColSums(
        *CEpetra::getVector(xID));
}

int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getRowMatrix(selfID)->RightScale(
        *CEpetra::getConstVector(xID));
}

boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->Filled();
}

double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormInf();
}

double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormOne();
}

int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalRows();
}

int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalCols();
}

int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyNonzeros();
}

int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRows();
}

int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyCols();
}

int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyDiagonals();
}

boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->LowerTriangular();
}

boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->UpperTriangular();
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstRowMatrix(
        selfID)->RowMatrixImporter());
}


} // extern "C"


//
// Definitions from CEpetra_RowMatrix_Cpp.hpp
//


/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
    return tableOfRowMatrixs().get(aid);
}

/* get Epetra_RowMatrix from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfRowMatrixs().get(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
    if (id.is_const) {
        return tableOfConstRowMatrixs().get(aid);
    } else {
        return tableOfRowMatrixs().get(aid);
    }
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstRowMatrixs().get(id);
    } else {
        return tableOfRowMatrixs().get(id);
    }
}

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
            tableOfRowMatrixs().storeShared(pobj));
}

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeConstRowMatrix( const Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
            tableOfConstRowMatrixs().storeShared(pobj));
}

/* dump contents of Epetra_RowMatrix and const Epetra_RowMatrix tables */
void
CEpetra::purgeRowMatrixTables(  )
{
    tableOfRowMatrixs().purge();
    tableOfConstRowMatrixs().purge();
}



