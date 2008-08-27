/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"

// ====================================================================== 
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix_in,
                            int OverlapLevel_in)  :
  Matrix_(Matrix_in),
  OverlapLevel_(OverlapLevel_in)
{
  // should not be here if no overlap
  if (OverlapLevel_in == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);
  
  NumMyRowsA_ = A().NumMyRows();

  // construct the external matrix
  vector<int> ExtElements; 

  Teuchos::RefCountPtr<Epetra_Map> TmpMap;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> TmpMatrix; 
  Teuchos::RefCountPtr<Epetra_Import> TmpImporter;

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 

  for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap) {
    if (TmpMatrix != Teuchos::null) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
    }

    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    vector<int> list(size); 

    int count = 0; 

    // define the set of rows that are in ColMap but not in RowMap 
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) { 
      int GID = ColMap->GID(i); 
      if (A().RowMatrixRowMap().LID(GID) == -1) { 
        vector<int>::iterator pos 
          = find(ExtElements.begin(),ExtElements.end(),GID); 
        if (pos == ExtElements.end()) { 
          ExtElements.push_back(GID);
          list[count] = GID; 
          ++count; 
        } 
      } 
    } 

    TmpMap = Teuchos::rcp( new Epetra_Map(-1,count, &list[0],0,Comm()) ); 

    TmpMatrix = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*TmpMap,0) ); 

    TmpImporter = Teuchos::rcp( new Epetra_Import(*TmpMap,A().RowMatrixRowMap()) ); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

  }

  // build the map containing all the nodes (original
  // matrix + extended matrix)
  vector<int> list(NumMyRowsA_ + ExtElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    list[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ExtElements.size() ; ++i)
    list[i + NumMyRowsA_] = ExtElements[i];

  Map_ = Teuchos::rcp( new Epetra_Map(-1, NumMyRowsA_ + ExtElements.size(),
				      &list[0], 0, Comm()) );
  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = Teuchos::rcp( new Epetra_Map(-1,ExtElements.size(),
					 &ExtElements[0],0,A().Comm()) );
  ExtMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*Map_,0) ); 

  ExtImporter_ = Teuchos::rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 
  ExtMatrix_->FillComplete(A().OperatorDomainMap(),*Map_);

  Importer_ = Teuchos::rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) );

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;
  NumMyCols_ = NumMyRows_;
  
  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
NumMyRowEntries(int MyRow, int & NumEntries) const
{
  if (MyRow < NumMyRowsA_)
    return(A().NumMyRowEntries(MyRow,NumEntries));
  else
    return(B().NumMyRowEntries(MyRow - NumMyRowsA_, NumEntries));
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, 
                 int * Indices) const
{
  int ierr;
  if (MyRow < NumMyRowsA_)
    ierr = A().ExtractMyRowCopy(MyRow,Length,NumEntries,Values,Indices);
  else
    ierr = B().ExtractMyRowCopy(MyRow - NumMyRowsA_,Length,NumEntries,
                                Values,Indices);

  IFPACK_RETURN(ierr);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(-1);
}


// ======================================================================
int Ifpack_OverlappingRowMatrix::
Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  vector<int> Ind(MaxNumEntries_);
  vector<double> Val(MaxNumEntries_);

  Y.PutScalar(0.0);

  // matvec with A (local rows)
  for (int i = 0 ; i < NumMyRowsA_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(A().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i] += Val[j] * X[k][Ind[j]];
      }
    }
  }

  // matvec with B (overlapping rows)
  for (int i = 0 ; i < NumMyRowsB_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(B().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i + NumMyRowsA_] += Val[j] * X[k][Ind[j]];
      }
    }
  }
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Multiply(UseTranspose(),X,Y));
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1);
}

// ======================================================================
Epetra_RowMatrix& Ifpack_OverlappingRowMatrix::B() const
{
  return(*ExtMatrix_);
}

// ======================================================================
const Epetra_BlockMap& Ifpack_OverlappingRowMatrix::Map() const
{
  return(*Map_);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ImportMultiVector(const Epetra_MultiVector& X, Epetra_MultiVector& OvX,
                  Epetra_CombineMode CM)
{
  OvX.Import(X,*Importer_,CM);
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExportMultiVector(const Epetra_MultiVector& OvX, Epetra_MultiVector& X,
                  Epetra_CombineMode CM)
{
  X.Export(OvX,*Importer_,CM);
  return(0);
}

