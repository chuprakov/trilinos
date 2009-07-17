#include "CTrilinos_config.h"

#include "Epetra_DataAccess.h"
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
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_CrsMatrix>& tableOfCrsMatrixs()
{
    static Table<Epetra_CrsMatrix>
        loc_tableOfCrsMatrixs(CT_Epetra_CrsMatrix_ID, "CT_Epetra_CrsMatrix_ID");
    return loc_tableOfCrsMatrixs;
}


} // namespace


//
// Definitions from CEpetra_CrsMatrix.h
//


extern "C" {


CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        CTrilinos::cast(tableOfCrsMatrixs(), id));
}

CTrilinos_Object_ID_t Epetra_CrsMatrix_Abstract ( 
  CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id);
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_Map> 
        pRowMap = CEpetra::getMap(RowMapID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *pRowMap, NumEntriesPerRow, StaticProfile)));

}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_Map> 
        pRowMap = CEpetra::getMap(RowMapID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *pRowMap, NumEntriesPerRow, StaticProfile)));

}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_Map> 
        pRowMap = CEpetra::getMap(RowMapID);
    const Teuchos::RCP<Epetra_Map> 
        pColMap = CEpetra::getMap(ColMapID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *pRowMap, *pColMap, NumEntriesPerRow, StaticProfile)));

}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<Epetra_Map> 
        pRowMap = CEpetra::getMap(RowMapID);
    const Teuchos::RCP<Epetra_Map> 
        pColMap = CEpetra::getMap(ColMapID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *pRowMap, *pColMap, NumEntriesPerRow, StaticProfile)));

}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  Epetra_DataAccess CV, CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<Epetra_CrsGraph> 
        pGraph = CEpetra::getCrsGraph(GraphID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *pGraph)));

}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<Epetra_CrsMatrix> 
        pMatrix = CEpetra::getCrsMatrix(MatrixID);

    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        *pMatrix)));

}

void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t id =
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(*selfID);
    tableOfCrsMatrixs().remove(&id);
    *selfID = CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(id);
}

void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID )
{
    Epetra_CrsMatrix& self = *( CEpetra::getCrsMatrix(selfID) );

    const Teuchos::RCP<Epetra_CrsMatrix> 
        psrc = CEpetra::getCrsMatrix(srcID);

    self = *psrc;
}

int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getCrsMatrix(selfID)->PutScalar(ScalarConstant);
}

int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getCrsMatrix(selfID)->Scale(ScalarConstant);
}

int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->InsertGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->SumIntoGlobalValues(
        GlobalRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->InsertMyValues(
        MyRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceMyValues(
        MyRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->SumIntoMyValues(
        MyRow, NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> 
        pDiagonal = CEpetra::getVector(DiagonalID);

    return CEpetra::getCrsMatrix(selfID)->ReplaceDiagonalValues(
        *pDiagonal);
}

int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage )
{
    return CEpetra::getCrsMatrix(selfID)->FillComplete(
        OptimizeDataStorage);
}

int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage )
{
    const Teuchos::RCP<Epetra_Map> 
        pDomainMap = CEpetra::getMap(DomainMapID);
    const Teuchos::RCP<Epetra_Map> 
        pRangeMap = CEpetra::getMap(RangeMapID);

    return CEpetra::getCrsMatrix(selfID)->FillComplete(
        *pDomainMap, *pRangeMap, OptimizeDataStorage);
}

int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->OptimizeStorage();
}

int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->MakeDataContiguous();
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractGlobalRowCopy(
        GlobalRow, Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractGlobalRowCopy(
        GlobalRow, Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> 
        pDiagonal = CEpetra::getVector(DiagonalID);

    return CEpetra::getCrsMatrix(selfID)->ExtractDiagonalCopy(
        *pDiagonal);
}

int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractGlobalRowView(
        GlobalRow, *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractMyRowView(
        MyRow, *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractGlobalRowView(
        GlobalRow, *NumEntries, *Values);
}

int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getCrsMatrix(selfID)->ExtractMyRowView(
        MyRow, *NumEntries, *Values);
}

int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);
    const Teuchos::RCP<Epetra_Vector> 
        py = CEpetra::getVector(yID);

    return CEpetra::getCrsMatrix(selfID)->Multiply(TransA, *px, *py);
}

int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);
    const Teuchos::RCP<Epetra_Vector> 
        py = CEpetra::getVector(yID);

    return CEpetra::getCrsMatrix(selfID)->Multiply1(TransA, *px, *py);
}

int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getCrsMatrix(selfID)->Multiply(TransA, *pX, *pY);
}

int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getCrsMatrix(selfID)->Multiply1(TransA, *pX, *pY);
}

int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);
    const Teuchos::RCP<Epetra_Vector> 
        py = CEpetra::getVector(yID);

    return CEpetra::getCrsMatrix(selfID)->Solve(
        Upper, Trans, UnitDiagonal, *px, *py);
}

int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getCrsMatrix(selfID)->Solve(
        Upper, Trans, UnitDiagonal, *pX, *pY);
}

int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->InvRowSums(*px);
}

int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->InvRowMaxs(*px);
}

int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->LeftScale(*px);
}

int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->InvColSums(*px);
}

int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->InvColMaxs(*px);
}

int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> 
        px = CEpetra::getVector(xID);

    return CEpetra::getCrsMatrix(selfID)->RightScale(*px);
}

boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->Filled();
}

boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->StorageOptimized();
}

boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->IndicesAreGlobal();
}

boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->IndicesAreLocal();
}

boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->IndicesAreContiguous();
}

boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->LowerTriangular();
}

boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->UpperTriangular();
}

boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NoDiagonal();
}

double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NormInf();
}

double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NormOne();
}

double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NormFrobenius();
}

int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumGlobalRows();
}

int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumGlobalCols();
}

int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyNonzeros();
}

int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyRows();
}

int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyCols();
}

int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyDiagonals();
}

int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getCrsMatrix(selfID)->NumGlobalEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getCrsMatrix(selfID)->NumAllocatedGlobalEntries(
        Row);
}

int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->MaxNumEntries();
}

int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->GlobalMaxNumEntries();
}

int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getCrsMatrix(selfID)->NumAllocatedMyEntries(Row);
}

int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->IndexBase();
}

boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->StaticGraph();
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeCrsGraph(
        &( CEpetra::getCrsMatrix(selfID)->Graph() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->RowMap() ));
}

int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pnewmap = CEpetra::getBlockMap(newmapID);

    return CEpetra::getCrsMatrix(selfID)->ReplaceRowMap(*pnewmap);
}

boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->HaveColMap();
}

int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pnewmap = CEpetra::getBlockMap(newmapID);

    return CEpetra::getCrsMatrix(selfID)->ReplaceColMap(*pnewmap);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->ColMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->DomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeImport(
        CEpetra::getCrsMatrix(selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeExport(
        CEpetra::getCrsMatrix(selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeComm(
        &( CEpetra::getCrsMatrix(selfID)->Comm() ));
}

int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CEpetra::getCrsMatrix(selfID)->LRID(GRID_in);
}

int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CEpetra::getCrsMatrix(selfID)->GRID(LRID_in);
}

int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CEpetra::getCrsMatrix(selfID)->LCID(GCID_in);
}

int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CEpetra::getCrsMatrix(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CEpetra::getCrsMatrix(selfID)->MyGRID(GRID_in);
}

boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CEpetra::getCrsMatrix(selfID)->MyLRID(LRID_in);
}

boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CEpetra::getCrsMatrix(selfID)->MyGCID(GCID_in);
}

boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CEpetra::getCrsMatrix(selfID)->MyLCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID )
{
    return CEpetra::getCrsMatrix(selfID)->MyGlobalRow(GID);
}

const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->Label();
}

int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in )
{
    return CEpetra::getCrsMatrix(selfID)->SetUseTranspose(
        UseTranspose_in);
}

int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getCrsMatrix(selfID)->Apply(*pX, *pY);
}

int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pX = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pY = CEpetra::getMultiVector(YID);

    return CEpetra::getCrsMatrix(selfID)->ApplyInverse(*pX, *pY);
}

boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->HasNormInf();
}

boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->UseTranspose();
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->OperatorRangeMap() ));
}

int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getCrsMatrix(selfID)->NumMyRowEntries(
        MyRow, *NumEntries);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeImport(
        CEpetra::getCrsMatrix(selfID)->RowMatrixImporter());
}

double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc )
{
    Epetra_CrsMatrix& self = *( CEpetra::getCrsMatrix(selfID) );

    return self[Loc];
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeMap(
        &( CEpetra::getCrsMatrix(selfID)->ImportMap() ));
}

int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->TransformToLocal();
}

int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID )
{
    const Teuchos::RCP<Epetra_Map> 
        pDomainMap = CEpetra::getMap(DomainMapID);
    const Teuchos::RCP<Epetra_Map> 
        pRangeMap = CEpetra::getMap(RangeMapID);

    return CEpetra::getCrsMatrix(selfID)->TransformToLocal(
        &(*pDomainMap), &(*pRangeMap));
}


} // extern "C"


//
// Definitions from CEpetra_CrsMatrix_Cpp.hpp
//


const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return tableOfCrsMatrixs().get(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfCrsMatrixs().get(id);
}

CT_Epetra_CrsMatrix_ID_t
CEpetra::storeCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().storeCopy(pobj));
}

void
CEpetra::purgeCrsMatrixTable(  )
{
    tableOfCrsMatrixs().purge();
}



