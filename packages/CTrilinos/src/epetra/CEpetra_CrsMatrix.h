#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_CRSMATRIX_H
#define CEPETRA_CRSMATRIX_H


#include "Epetra_DataAccess.h"
#include "CEpetra_Map.h"
#include "CEpetra_CrsGraph.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Import.h"
#include "CEpetra_Export.h"
#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_CrsMatrix_Abstract ( 
  CT_Epetra_CrsMatrix_ID_t id );

/* Original C++ prototype:
   Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  Epetra_DataAccess CV, CT_Epetra_CrsGraph_ID_t GraphID );

/* Original C++ prototype:
   Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix);
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID );

/* Original C++ prototype:
   virtual ~Epetra_CrsMatrix();
*/
void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID );

/* Original C++ prototype:
   Epetra_CrsMatrix& operator=(const Epetra_CrsMatrix& src);
*/
void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID );

/* Original C++ prototype:
   int PutScalar(double ScalarConstant);
*/
int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

/* Original C++ prototype:
   int Scale(double ScalarConstant);
*/
int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

/* Original C++ prototype:
   virtual int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   virtual int ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   virtual int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   int InsertMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   int ReplaceMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   int SumIntoMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
*/
int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/* Original C++ prototype:
   int ReplaceDiagonalValues(const Epetra_Vector& Diagonal);
*/
int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/* Original C++ prototype:
   int FillComplete(bool OptimizeDataStorage = true);
*/
int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage );

/* Original C++ prototype:
   int FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap, bool OptimizeDataStorage = true);
*/
int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage );

/* Original C++ prototype:
   int OptimizeStorage();
*/
int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int MakeDataContiguous();
*/
int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, int* Indices) const;
*/
int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/* Original C++ prototype:
   int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, int* Indices) const;
*/
int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/* Original C++ prototype:
   int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const;
*/
int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values );

/* Original C++ prototype:
   int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values) const;
*/
int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values );

/* Original C++ prototype:
   int ExtractDiagonalCopy(Epetra_Vector& Diagonal) const;
*/
int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/* Original C++ prototype:
   int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values, int*& Indices) const;
*/
int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices );

/* Original C++ prototype:
   int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const;
*/
int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices );

/* Original C++ prototype:
   int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values) const;
*/
int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values );

/* Original C++ prototype:
   int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values) const;
*/
int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values );

/* Original C++ prototype:
   int Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;
*/
int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

/* Original C++ prototype:
   int Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;
*/
int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

/* Original C++ prototype:
   int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   int Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;
*/
int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID );

/* Original C++ prototype:
   int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   int InvRowSums(Epetra_Vector& x) const;
*/
int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   int InvRowMaxs(Epetra_Vector& x) const;
*/
int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   int LeftScale(const Epetra_Vector& x);
*/
int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   int InvColSums(Epetra_Vector& x) const;
*/
int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   int InvColMaxs(Epetra_Vector& x) const;
*/
int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   int RightScale(const Epetra_Vector& x);
*/
int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   bool Filled() const;
*/
boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool StorageOptimized() const;
*/
boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool IndicesAreGlobal() const;
*/
boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool IndicesAreLocal() const;
*/
boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool IndicesAreContiguous() const;
*/
boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool LowerTriangular() const;
*/
boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool UpperTriangular() const;
*/
boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool NoDiagonal() const;
*/
boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   double NormInf() const;
*/
double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   double NormOne() const;
*/
double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   double NormFrobenius() const;
*/
double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalNonzeros() const;
*/
int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalRows() const;
*/
int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalCols() const;
*/
int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalDiagonals() const;
*/
int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyNonzeros() const;
*/
int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyRows() const;
*/
int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyCols() const;
*/
int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyDiagonals() const;
*/
int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalEntries(int Row) const;
*/
int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/* Original C++ prototype:
   int NumAllocatedGlobalEntries(int Row) const;
*/
int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/* Original C++ prototype:
   int MaxNumEntries() const;
*/
int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int GlobalMaxNumEntries() const;
*/
int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyEntries(int Row) const;
*/
int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/* Original C++ prototype:
   int NumAllocatedMyEntries(int Row) const;
*/
int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/* Original C++ prototype:
   int IndexBase() const;
*/
int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool StaticGraph();
*/
boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_CrsGraph& Graph() const;
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& RowMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int ReplaceRowMap(const Epetra_BlockMap& newmap);
*/
int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/* Original C++ prototype:
   bool HaveColMap() const;
*/
boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int ReplaceColMap(const Epetra_BlockMap& newmap);
*/
int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/* Original C++ prototype:
   const Epetra_Map& ColMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& DomainMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& RangeMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Import* Importer() const;
*/
CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Export* Exporter() const;
*/
CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Comm& Comm() const;
*/
CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int LRID( int GRID_in) const;
*/
int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

/* Original C++ prototype:
   int GRID( int LRID_in) const;
*/
int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

/* Original C++ prototype:
   int LCID( int GCID_in) const;
*/
int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

/* Original C++ prototype:
   int GCID( int LCID_in) const;
*/
int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

/* Original C++ prototype:
   bool MyGRID(int GRID_in) const;
*/
boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

/* Original C++ prototype:
   bool MyLRID(int LRID_in) const;
*/
boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

/* Original C++ prototype:
   bool MyGCID(int GCID_in) const;
*/
boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

/* Original C++ prototype:
   bool MyLCID(int LCID_in) const;
*/
boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

/* Original C++ prototype:
   bool MyGlobalRow(int GID) const;
*/
boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID );

/* Original C++ prototype:
   const char* Label() const;
*/
const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int SetUseTranspose(bool UseTranspose_in);
*/
int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in );

/* Original C++ prototype:
   int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   bool HasNormInf() const;
*/
boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   bool UseTranspose() const;
*/
boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& OperatorDomainMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& OperatorRangeMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int NumMyRowEntries(int MyRow, int& NumEntries) const;
*/
int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries );

/* Original C++ prototype:
   const Epetra_Map& RowMatrixRowMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Map& RowMatrixColMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Import* RowMatrixImporter() const;
*/
CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   inline double* operator[] (int Loc) const;
*/
double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc );

/* Original C++ prototype:
   const Epetra_Map& ImportMap() const;
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int TransformToLocal();
*/
int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/* Original C++ prototype:
   int TransformToLocal(const Epetra_Map* DomainMap, const Epetra_Map* RangeMap);
*/
int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_CRSMATRIX_H */

