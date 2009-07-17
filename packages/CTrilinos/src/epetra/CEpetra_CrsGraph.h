#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_CRSGRAPH_H
#define CEPETRA_CRSGRAPH_H


#include "Epetra_DataAccess.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Import.h"
#include "CEpetra_Export.h"
#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_CrsGraph_Abstract ( 
  CT_Epetra_CrsGraph_ID_t id );

/* Original C++ prototype:
   Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const int* NumIndicesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, const int* NumIndicesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, int NumIndicesPerRow, bool StaticProfile = false);
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile );

/* Original C++ prototype:
   Epetra_CrsGraph(const Epetra_CrsGraph& Graph);
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID );

/* Original C++ prototype:
   virtual ~Epetra_CrsGraph();
*/
void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID );

/* Original C++ prototype:
   int InsertGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
*/
int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );

/* Original C++ prototype:
   int RemoveGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
*/
int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );

/* Original C++ prototype:
   int RemoveGlobalIndices(int Row);
*/
int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int InsertMyIndices(int LocalRow, int NumIndices, int* Indices);
*/
int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );

/* Original C++ prototype:
   int RemoveMyIndices(int LocalRow, int NumIndices, int* Indices);
*/
int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );

/* Original C++ prototype:
   int RemoveMyIndices(int Row);
*/
int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int FillComplete();
*/
int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap);
*/
int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );

/* Original C++ prototype:
   int OptimizeStorage();
*/
int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices, int& NumIndices, int* Indices) const;
*/
int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );

/* Original C++ prototype:
   int ExtractMyRowCopy(int LocalRow, int LenOfIndices, int& NumIndices, int* Indices) const;
*/
int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );

/* Original C++ prototype:
   int ExtractGlobalRowView(int GlobalRow, int& NumIndices, int*& Indices) const;
*/
int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices );

/* Original C++ prototype:
   int ExtractMyRowView(int LocalRow, int& NumIndices, int*& Indices) const;
*/
int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices );

/* Original C++ prototype:
   bool Filled() const;
*/
boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool StorageOptimized() const;
*/
boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool IndicesAreGlobal() const;
*/
boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool IndicesAreLocal() const;
*/
boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool LowerTriangular() const;
*/
boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool UpperTriangular() const;
*/
boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool NoDiagonal() const;
*/
boolean Epetra_CrsGraph_NoDiagonal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   bool MyGlobalRow(int GID) const;
*/
boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID );

/* Original C++ prototype:
   bool HaveColMap() const;
*/
boolean Epetra_CrsGraph_HaveColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyRows() const;
*/
int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalRows() const;
*/
int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyCols() const;
*/
int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalCols() const;
*/
int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalNonzeros() const;
*/
int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalDiagonals() const;
*/
int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyDiagonals() const;
*/
int Epetra_CrsGraph_NumMyDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyBlockRows() const;
*/
int Epetra_CrsGraph_NumMyBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalBlockRows() const;
*/
int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyBlockCols() const;
*/
int Epetra_CrsGraph_NumMyBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalBlockCols() const;
*/
int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyBlockDiagonals() const;
*/
int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalBlockDiagonals() const;
*/
int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalEntries() const;
*/
int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyEntries() const;
*/
int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int MaxRowDim() const;
*/
int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int GlobalMaxRowDim() const;
*/
int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int MaxColDim() const;
*/
int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int GlobalMaxColDim() const;
*/
int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyNonzeros() const;
*/
int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumGlobalIndices(int Row) const;
*/
int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int NumAllocatedGlobalIndices(int Row) const;
*/
int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int MaxNumIndices() const;
*/
int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int GlobalMaxNumIndices() const;
*/
int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int MaxNumNonzeros() const;
*/
int Epetra_CrsGraph_MaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int GlobalMaxNumNonzeros() const;
*/
int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int NumMyIndices(int Row) const;
*/
int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int NumAllocatedMyIndices(int Row) const;
*/
int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/* Original C++ prototype:
   int IndexBase() const;
*/
int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_BlockMap& RowMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int ReplaceRowMap(const Epetra_BlockMap& newmap);
*/
int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/* Original C++ prototype:
   int ReplaceColMap(const Epetra_BlockMap& newmap);
*/
int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/* Original C++ prototype:
   const Epetra_BlockMap& ColMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_BlockMap& DomainMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_BlockMap& RangeMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Import* Importer() const;
*/
CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Export* Exporter() const;
*/
CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   const Epetra_Comm& Comm() const;
*/
CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int LRID(int GRID_in) const;
*/
int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

/* Original C++ prototype:
   int GRID(int LRID_in) const;
*/
int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

/* Original C++ prototype:
   int LCID(int GCID_in) const;
*/
int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

/* Original C++ prototype:
   int GCID(int LCID_in) const;
*/
int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

/* Original C++ prototype:
   bool MyGRID(int GRID_in) const;
*/
boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

/* Original C++ prototype:
   bool MyLRID(int LRID_in) const;
*/
boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

/* Original C++ prototype:
   bool MyGCID(int GCID_in) const;
*/
boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

/* Original C++ prototype:
   bool MyLCID(int LCID_in) const;
*/
boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

/* Original C++ prototype:
   inline int* operator[]( int Loc ) const;
*/
int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc );

/* Original C++ prototype:
   Epetra_CrsGraph& operator = (const Epetra_CrsGraph& Source);
*/
void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID );

/* Original C++ prototype:
   const Epetra_BlockMap& ImportMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int TransformToLocal();
*/
int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/* Original C++ prototype:
   int TransformToLocal(const Epetra_BlockMap* DomainMap, const Epetra_BlockMap* RangeMap);
*/
int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_CRSGRAPH_H */

