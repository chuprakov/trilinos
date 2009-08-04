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
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_CrsMatrix */
Table<Epetra_CrsMatrix>& tableOfCrsMatrixs()
{
    static Table<Epetra_CrsMatrix>
        loc_tableOfCrsMatrixs(CT_Epetra_CrsMatrix_ID, "CT_Epetra_CrsMatrix_ID", false);
    return loc_tableOfCrsMatrixs;
}

/* table to hold objects of type const Epetra_CrsMatrix */
Table<const Epetra_CrsMatrix>& tableOfConstCrsMatrixs()
{
    static Table<const Epetra_CrsMatrix>
        loc_tableOfConstCrsMatrixs(CT_Epetra_CrsMatrix_ID, "CT_Epetra_CrsMatrix_ID", true);
    return loc_tableOfConstCrsMatrixs;
}


} // namespace


//
// Definitions from CEpetra_CrsMatrix.h
//


extern "C" {


CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstCrsMatrixs(), id);
    } else {
        newid = CTrilinos::cast(tableOfCrsMatrixs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(newid);
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
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *CEpetra::getConstMap(RowMapID), NumEntriesPerRow, 
        StaticProfile)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *CEpetra::getConstMap(RowMapID), NumEntriesPerRow, 
        StaticProfile)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *CEpetra::getConstMap(RowMapID), 
        *CEpetra::getConstMap(ColMapID), NumEntriesPerRow, 
        StaticProfile)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *CEpetra::getConstMap(RowMapID), 
        *CEpetra::getConstMap(ColMapID), NumEntriesPerRow, 
        StaticProfile)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  Epetra_DataAccess CV, CT_Epetra_CrsGraph_ID_t GraphID )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        CV, *CEpetra::getConstCrsGraph(GraphID))));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store(new Epetra_CrsMatrix(
        *CEpetra::getConstCrsMatrix(MatrixID))));
}

void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstCrsMatrixs().remove(&aid);
    } else {
        tableOfCrsMatrixs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(aid);
}

void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID )
{
    Epetra_CrsMatrix& self = *( CEpetra::getCrsMatrix(selfID) );

    self = *CEpetra::getConstCrsMatrix(srcID);
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
    return CEpetra::getCrsMatrix(selfID)->ReplaceDiagonalValues(
        *CEpetra::getConstVector(DiagonalID));
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
    return CEpetra::getCrsMatrix(selfID)->FillComplete(
        *CEpetra::getConstMap(DomainMapID), 
        *CEpetra::getConstMap(RangeMapID), OptimizeDataStorage);
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
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowCopy(
        GlobalRow, Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowCopy(
        GlobalRow, Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowCopy(
        MyRow, Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractDiagonalCopy(
        *CEpetra::getVector(DiagonalID));
}

int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowView(
        GlobalRow, *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowView(
        MyRow, *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowView(
        GlobalRow, *NumEntries, *Values);
}

int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowView(
        MyRow, *NumEntries, *Values);
}

int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Multiply(
        TransA, *CEpetra::getConstVector(xID), 
        *CEpetra::getVector(yID));
}

int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Multiply1(
        TransA, *CEpetra::getConstVector(xID), 
        *CEpetra::getVector(yID));
}

int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Multiply(
        TransA, *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Multiply1(
        TransA, *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Solve(
        Upper, Trans, UnitDiagonal, *CEpetra::getConstVector(xID), 
        *CEpetra::getVector(yID));
}

int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Solve(
        Upper, Trans, UnitDiagonal, 
        *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstCrsMatrix(selfID)->InvRowSums(
        *CEpetra::getVector(xID));
}

int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstCrsMatrix(selfID)->InvRowMaxs(
        *CEpetra::getVector(xID));
}

int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getCrsMatrix(selfID)->LeftScale(
        *CEpetra::getConstVector(xID));
}

int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstCrsMatrix(selfID)->InvColSums(
        *CEpetra::getVector(xID));
}

int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getConstCrsMatrix(selfID)->InvColMaxs(
        *CEpetra::getVector(xID));
}

int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    return CEpetra::getCrsMatrix(selfID)->RightScale(
        *CEpetra::getConstVector(xID));
}

boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Filled();
}

boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->StorageOptimized();
}

boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->IndicesAreGlobal();
}

boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->IndicesAreLocal();
}

boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->IndicesAreContiguous();
}

boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->LowerTriangular();
}

boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->UpperTriangular();
}

boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NoDiagonal();
}

double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormInf();
}

double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormOne();
}

double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormFrobenius();
}

int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalRows();
}

int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalCols();
}

int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyNonzeros();
}

int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyRows();
}

int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyCols();
}

int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyDiagonals();
}

int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumAllocatedGlobalEntries(
        Row);
}

int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->MaxNumEntries();
}

int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->GlobalMaxNumEntries();
}

int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumAllocatedMyEntries(
        Row);
}

int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->IndexBase();
}

boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->StaticGraph();
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstCrsGraph(
        &( CEpetra::getConstCrsMatrix(selfID)->Graph() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->RowMap() ));
}

int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceRowMap(
        *CEpetra::getConstBlockMap(newmapID));
}

boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->HaveColMap();
}

int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceColMap(
        *CEpetra::getConstBlockMap(newmapID));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->ColMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->DomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(
        CEpetra::getConstCrsMatrix(selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstExport(
        CEpetra::getConstCrsMatrix(selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstComm(
        &( CEpetra::getConstCrsMatrix(selfID)->Comm() ));
}

int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->LRID(GRID_in);
}

int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->GRID(LRID_in);
}

int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->LCID(GCID_in);
}

int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->MyGRID(GRID_in);
}

boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->MyLRID(LRID_in);
}

boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->MyGCID(GCID_in);
}

boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->MyLCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID )
{
    return CEpetra::getConstCrsMatrix(selfID)->MyGlobalRow(GID);
}

const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Label();
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
    return CEpetra::getConstCrsMatrix(selfID)->Apply(
        *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    return CEpetra::getConstCrsMatrix(selfID)->ApplyInverse(
        *CEpetra::getConstMultiVector(XID), 
        *CEpetra::getMultiVector(YID));
}

boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->HasNormInf();
}

boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->UseTranspose();
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->OperatorRangeMap() ));
}

int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyRowEntries(
        MyRow, *NumEntries);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(
        CEpetra::getConstCrsMatrix(selfID)->RowMatrixImporter());
}

double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc )
{
    const Epetra_CrsMatrix& self = *( CEpetra::getConstCrsMatrix(selfID) );

    return self[Loc];
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(
        &( CEpetra::getConstCrsMatrix(selfID)->ImportMap() ));
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
    return CEpetra::getCrsMatrix(selfID)->TransformToLocal(
        CEpetra::getConstMap(DomainMapID).getRawPtr(), 
        CEpetra::getConstMap(RangeMapID).getRawPtr());
}


} // extern "C"


//
// Definitions from CEpetra_CrsMatrix_Cpp.hpp
//


/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id);
    return tableOfCrsMatrixs().get(aid);
}

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CTrilinos_Object_ID_t id )
{
    return tableOfCrsMatrixs().get(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id);
    if (id.is_const) {
        return tableOfConstCrsMatrixs().get(aid);
    } else {
        return tableOfCrsMatrixs().get(aid);
    }
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstCrsMatrixs().get(id);
    } else {
        return tableOfCrsMatrixs().get(id);
    }
}

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
            tableOfCrsMatrixs().storeCopy(pobj));
}

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeConstCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
            tableOfConstCrsMatrixs().storeCopy(pobj));
}

/* dump contents of Epetra_CrsMatrix and const Epetra_CrsMatrix tables */
void
CEpetra::purgeCrsMatrixTables(  )
{
    tableOfCrsMatrixs().purge();
    tableOfConstCrsMatrixs().purge();
}



