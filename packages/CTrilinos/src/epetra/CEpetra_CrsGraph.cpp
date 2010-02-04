
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

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_CrsGraph.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_CrsGraph.h
//


extern "C" {


CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CrsGraph_Generalize ( 
  CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RowMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(new Epetra_CrsGraph((Epetra_DataAccess) CV, 
        *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RowMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(new Epetra_CrsGraph((Epetra_DataAccess) CV, 
        *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(ColMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(new Epetra_CrsGraph((Epetra_DataAccess) CV, 
        *RowMap, *ColMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(ColMapID);
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(new Epetra_CrsGraph((Epetra_DataAccess) CV, 
        *RowMap, *ColMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = 
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(GraphID);
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(new Epetra_CrsGraph(*Graph));
}

void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->InsertGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RemoveGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RemoveGlobalIndices(Row);
}

int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->InsertMyIndices(LocalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RemoveMyIndices(LocalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RemoveMyIndices(Row);
}

int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->FillComplete();
}

int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RangeMapID);
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->FillComplete(*DomainMap, *RangeMap);
}

int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->OptimizeStorage();
}

int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ExtractMyRowCopy(LocalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumIndices, *Indices);
}

int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ExtractMyRowView(LocalRow, 
        *NumIndices, *Indices);
}

boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->Filled()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->StorageOptimized()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->IndicesAreGlobal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->IndicesAreLocal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->UpperTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_NoDiagonal ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NoDiagonal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MyGlobalRow(GID)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_HaveColMap ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->HaveColMap()) ? TRUE : FALSE);
}

int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyRows();
}

int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalRows();
}

int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyCols();
}

int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalCols();
}

int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsGraph_NumMyDiagonals ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyDiagonals();
}

int Epetra_CrsGraph_NumMyBlockRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyBlockRows();
}

int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalBlockRows();
}

int Epetra_CrsGraph_NumMyBlockCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyBlockCols();
}

int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalBlockCols();
}

int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalEntries();
}

int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyEntries();
}

int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MaxRowDim();
}

int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GlobalMaxRowDim();
}

int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MaxColDim();
}

int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GlobalMaxColDim();
}

int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyNonzeros();
}

int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumGlobalIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumAllocatedGlobalIndices(Row);
}

int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MaxNumIndices();
}

int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GlobalMaxNumIndices();
}

int Epetra_CrsGraph_MaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MaxNumNonzeros();
}

int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GlobalMaxNumNonzeros();
}

int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumMyIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->NumAllocatedMyIndices(Row);
}

int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->IndexBase();
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RowMap() ));
}

int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(newmapID);
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ReplaceRowMap(*newmap);
}

int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(newmapID);
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ReplaceColMap(*newmap);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ColMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->DomainMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Import, CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->Comm() ));
}

int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->LRID(GRID_in);
}

int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GRID(LRID_in);
}

int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->LCID(GCID_in);
}

int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MyGRID(GRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MyLRID(LRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MyGCID(GCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->MyLCID(LCID_in)) ? TRUE : FALSE);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->ImportMap() ));
}

int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->TransformToLocal();
}

int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(RangeMapID);
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID)->TransformToLocal(
        DomainMap.getRawPtr(), RangeMap.getRawPtr());
}

int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc )
{
    const Epetra_CrsGraph& self = *( 
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID) );

    return self[Loc];
}

void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID )
{
    Epetra_CrsGraph& self = *( CTrilinos::tableRepos().get<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(selfID) );

    const Teuchos::RCP<const Epetra_CrsGraph> Source = 
        CTrilinos::tableRepos().getConst<Epetra_CrsGraph, 
        CT_Epetra_CrsGraph_ID_t>(SourceID);
    self = *Source;
}


} // extern "C"


//
// Definitions from CEpetra_CrsGraph_Cpp.hpp
//


/* get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, CT_Epetra_CrsGraph_ID_t>(id);
}

/* get Epetra_CrsGraph from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, CT_Epetra_CrsGraph_ID_t>(id);
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeCrsGraph( Epetra_CrsGraph *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, CT_Epetra_CrsGraph_ID_t>(pobj, false);
}

/* store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeConstCrsGraph( const Epetra_CrsGraph *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_CrsGraph, CT_Epetra_CrsGraph_ID_t>(pobj, false);
}



