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
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_RowMatrix>& tableOfRowMatrixs()
{
    static Table<Epetra_RowMatrix> loc_tableOfRowMatrixs(CT_Epetra_RowMatrix_ID, "CT_Epetra_RowMatrix_ID");
    return loc_tableOfRowMatrixs;
}


} // namespace


//
// Definitions from CEpetra_RowMatrix.h
//


extern "C" {


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfRowMatrixs(), id);
}

void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID )
{
    tableOfRowMatrixs().remove(selfID);
}

int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getRowMatrix(selfID)->NumMyRowEntries(
        MyRow, *NumEntries);
}

int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->MaxNumEntries();
}

int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getRowMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values, Indices);
}

int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> 
        pDiagonal = CEpetra::getVector(DiagonalID);

    return CEpetra::getRowMatrix(selfID)->ExtractDiagonalCopy(
        *pDiagonal);
}

int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getRowMatrix(selfID)->Multiply(TransA, *pX, *pY);
}

int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getRowMatrix(selfID)->Solve(
        Upper, Trans, UnitDiagonal, *pX, *pY);
}

int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getRowMatrix(selfID)->InvRowSums(*px);
}

int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getRowMatrix(selfID)->LeftScale(*px);
}

int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getRowMatrix(selfID)->InvColSums(*px);
}

int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getRowMatrix(selfID)->RightScale(*px);
}

boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->Filled();
}

double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NormInf();
}

double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NormOne();
}

int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumGlobalRows();
}

int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumGlobalCols();
}

int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumMyNonzeros();
}

int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumMyRows();
}

int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumMyCols();
}

int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->NumMyDiagonals();
}

boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->LowerTriangular();
}

boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getRowMatrix(selfID)->UpperTriangular();
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getRowMatrix(selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getRowMatrix(selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeImport(
        CEpetra::getRowMatrix(selfID)->RowMatrixImporter());
}


} // extern "C"


//
// Definitions from CEpetra_RowMatrix_Cpp.hpp
//


const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return tableOfRowMatrixs().get(id);
}

CT_Epetra_RowMatrix_ID_t
CEpetra::storeRowMatrix( const Epetra_RowMatrix *pobj )
{
    return tableOfRowMatrixs().storeCopy(pobj);
}

void
CEpetra::purgeRowMatrixTable(  )
{
    tableOfRowMatrixs().purge();
}



