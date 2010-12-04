//@HEADER
// ************************************************************************
// 
//               EpetraExt: Extended Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const int * RowIndices,
	int NumBlockRows,
        const Epetra_Comm & GlobalComm ) 
{
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map *GlobalMap = 
    new Epetra_Map( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, 
		    GlobalComm );

  return GlobalMap;
}

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const std::vector<int>& RowIndices,
	const Epetra_Comm & GlobalComm ) 
{
  return GenerateBlockMap(BaseMap, &RowIndices[0], RowIndices.size(), 
			  GlobalComm);
}

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const Epetra_BlockMap& BlockMap,
	const Epetra_Comm & GlobalComm ) 
{
  return GenerateBlockMap(BaseMap, 
			  BlockMap.MyGlobalElements(), 
			  BlockMap.NumMyElements(), 
			  GlobalComm);
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm ) 
{

  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int BaseRow = BaseMap.GID(j);
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < StencilSize; ++k )
      {
        int ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm ) 
{

  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices(MaxIndices);
  vector<double> Values(MaxIndices); 
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseMatrix.ExtractMyRowCopy( j, MaxIndices, NumIndices, &Values[0], &Indices[0] );
      for( int l = 0; l < NumIndices; ++l ) Indices[l] = BaseColMap.GID(Indices[l]);

      for( int k = 0; k < StencilSize; ++k )
      {
        int ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm ) 
{
  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  const Epetra_BlockMap & BlockRowMap = LocalBlockGraph.RowMap();
  const Epetra_BlockMap & BlockColMap = LocalBlockGraph.ColMap();
  int NumBlockRows = BlockRowMap.NumMyElements();
  vector<int> RowIndices(NumBlockRows);
  BlockRowMap.MyGlobalElements(&RowIndices[0]);
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  int NumBlockIndices, NumBaseIndices;
  int *BlockIndices, *BaseIndices;
  for( int i = 0; i < NumBlockRows; ++i )
  {
    LocalBlockGraph.ExtractMyRowView(i, NumBlockIndices, BlockIndices);
    
    for( int j = 0; j < Size; ++j )
    {
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseGraph.ExtractMyRowView( j, NumBaseIndices, BaseIndices );
      for( int k = 0; k < NumBlockIndices; ++k )
      {
        int ColOffset = BlockColMap.GID(BlockIndices[k]) * Offset;

        for( int l = 0; l < NumBaseIndices; ++l )
          Indices[l] = BaseGraph.GCID(BaseIndices[l]) + ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumBaseIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

//==============================================================================
void BlockUtility::GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, 
				      std::vector<int> RowIndices, 
				      std::vector< std::vector<int> >& RowStencil)
{
  // Get row indices
  int NumMyRows = LocalBlockGraph.NumMyRows();
  RowIndices.resize(NumMyRows);
  const Epetra_BlockMap& RowMap = LocalBlockGraph.RowMap();
  RowMap.MyGlobalElements(&RowIndices[0]);

  // Get stencil
  RowStencil.resize(NumMyRows);
  if (LocalBlockGraph.IndicesAreGlobal()) {
    for (int i=0; i<NumMyRows; i++) {
      int Row = RowIndices[i];
      int NumCols = LocalBlockGraph.NumGlobalIndices(Row);
      RowStencil[i].resize(NumCols);
      LocalBlockGraph.ExtractGlobalRowCopy(Row, NumCols, NumCols, 
					   &RowStencil[i][0]);
      for (int k=0; k<NumCols; k++)
	RowStencil[i][k] -= Row;
    }
  }
  else {
    for (int i=0; i<NumMyRows; i++) {
      int NumCols = LocalBlockGraph.NumMyIndices(i);
      RowStencil[i].resize(NumCols);
      LocalBlockGraph.ExtractMyRowCopy(i, NumCols, NumCols, 
				       &RowStencil[i][0]);
      for (int k=0; k<NumCols; k++)
	RowStencil[i][k] = LocalBlockGraph.GCID(RowStencil[i][k]) - 
	  RowIndices[i];
    }
  }
}

//==============================================================================
int BlockUtility::CalculateOffset(const Epetra_BlockMap & BaseMap)
{
  int MaxGID = BaseMap.MaxAllGID();

//   int Offset = 1;
//   while( Offset <= MaxGID ) Offset *= 10;

//   return Offset;

  return MaxGID+1;
}


} //namespace EpetraExt
