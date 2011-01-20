//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_MATRIXMATRIX_DEF_HPP
#define TPETRA_MATRIXMATRIX_DEF_HPP

#include "Tpetra_MatrixMatrix_decl.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_map.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MMHelpers_def.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include <algorithm>

#ifdef DOXYGEN_USE_ONLY
  //#include "Tpetra_MMMultiply_decl.hpp"
#endif

/*! \file Tpetra_MMMultiply_def.hpp 

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */

namespace Tpetra {

//kernel method for computing the local portion of C = A*B
template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::mult_A_B(
  RCP<CrsMatrixStruct_t >& Aview, 
  RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{
  LocalOrdinal C_firstCol = Bview->colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview->colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = OrdinalTraits<LocalOrdinal>::invalid();

  //int* bcols = Bview.colMap->MyGlobalElements();
  ArrayView<const GlobalOrdinal> bcols =Bview->colMap->getNodeElementList();
  //int* bcols_import = NULL;
  ArrayView<const GlobalOrdinal> bcols_import = null;
  if (Bview->importColMap != null) {
    C_firstCol_import = Bview->importColMap->getMinLocalIndex();
    C_lastCol_import = Bview->importColMap->getMaxLocalIndex();

    bcols_import = Bview->importColMap->getNodeElementList();
  }

  size_t C_numCols = C_lastCol - C_firstCol + OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;
  Array<Scalar> dwork = Array<Scalar>(C_numCols);
  //int* iwork = new int[C_numCols];
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);

  Array<Scalar> C_row_i = dwork;
  //int* C_cols = iwork;
  Array<GlobalOrdinal> C_cols = iwork;

  //int C_row_i_length, i, j, k;
  size_t C_row_i_length, i, j, k;

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.

  bool C_filled = C.isFillComplete();

  //loop over the rows of A.
  for(i=0; i<Aview->numRows; ++i) {

    //only navigate the local portion of Aview... (It's probable that we
    //imported more of A than we need for A*B, because other cases like A^T*B 
    //need the extra rows.)
    if (Aview->remote[i]) {
      continue;
    }

    ArrayView<const LocalOrdinal> Aindices_i = Aview->indices[i];
    ArrayView<const Scalar> Aval_i  = Aview->values[i];

    GlobalOrdinal global_row = Aview->rowMap->getGlobalElement(i);


    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.

    for(k=OrdinalTraits<size_t>::zero(); k<Aview->numEntriesPerRow[i]; ++k) {
      LocalOrdinal Ak = Bview->rowMap->getLocalElement(Aview->colMap->getGlobalElement(Aindices_i[k]));
      Scalar Aval = Aval_i[k];

      ArrayView<const LocalOrdinal> Bcol_inds = Bview->indices[Ak];
      ArrayView<const Scalar> Bvals_k = Bview->values[Ak];

      C_row_i_length = OrdinalTraits<size_t>::zero();

      if (Bview->remote[Ak]) {
        for(j=OrdinalTraits<size_t>::zero(); j<Bview->numEntriesPerRow[Ak]; ++j) {
          C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          C_cols[C_row_i_length++] = bcols_import[Bcol_inds[j]];
        }
      }
      else {
        for(j=OrdinalTraits<size_t>::zero(); j<Bview->numEntriesPerRow[Ak]; ++j) {
          C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          C_cols[C_row_i_length++] = bcols[Bcol_inds[j]];
        }
      }
	  /*
	  std::cout << "About to insert row: " << global_row << std::endl;
	  ArrayView<const Scalar> C_row_iView = C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length);
      typename ArrayView<const Scalar>::iterator it = C_row_iView.begin();
      for(; it != C_row_iView.end(); ++it){
        std::cout << *it << ", ";
      }
      std::cout << std::endl;*/

      //
      //Now put the C_row_i values into C.
      //

      C_filled ?
        C.sumIntoGlobalValues(global_row, C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length), C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length))
        :
        C.insertGlobalValues(global_row, C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length), C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length));

    }
  }

  //delete [] dwork;
  //delete [] iwork;

  return(0);
}

//kernel method for computing the local portion of C = A*B^T
template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::mult_A_Btrans(
  const RCP<const CrsMatrixStruct_t> & Aview, 
  const RCP<const CrsMatrixStruct_t> & Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node> & C)
{
  size_t maxlen = 0;
  for (size_t i=0; i<Aview->numRows; ++i) {
    if (Aview->numEntriesPerRow[i] > maxlen) maxlen = Aview->numEntriesPerRow[i];
  }
  for (size_t i=0; i<Bview->numRows; ++i) {
    if (Bview->numEntriesPerRow[i] > maxlen) maxlen = Bview->numEntriesPerRow[i];
  }

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);

  size_t numBrows = Bview->numRows;

  Array<GlobalOrdinal>     iwork(maxlen*2);
  ArrayView<GlobalOrdinal>  Aind = iwork(0,maxlen);
  ArrayView<GlobalOrdinal>  Bind = iwork(maxlen,maxlen);

  ArrayView<const GlobalOrdinal> bgids = Bview->colMap->getNodeElementList();

  Array<Scalar>      vals(maxlen*2);
  ArrayView<Scalar> bvals = vals(0,maxlen);
  ArrayView<Scalar> avals = vals(maxlen, maxlen);

  const GlobalOrdinal max_all_b = Bview->colMap->getMaxAllGlobalIndex();
  const GlobalOrdinal min_all_b = Bview->colMap->getMinAllGlobalIndex();

  // next create arrays indicating the first and last column-index in
  // each row of B, so that we can know when to skip certain rows below.
  // This will provide a large performance gain for banded matrices, and
  // a somewhat smaller gain for *most* other matrices.
  Array<GlobalOrdinal> b_firstcol(numBrows);
  Array<GlobalOrdinal> b_lastcol(numBrows);
  GlobalOrdinal temp;
  for(size_t i=0; i<numBrows; ++i) {
    b_firstcol[i] = max_all_b;
    b_lastcol[i] = min_all_b;

    size_t Blen_i = Bview->numEntriesPerRow[i];
    if (Blen_i < 1) continue;
    ArrayView<const LocalOrdinal> Bindices_i = Bview->indices[i];

    if (Bview->remote[i]) {
      for(size_t k=0; k<Blen_i; ++k) {
        temp = Bview->importColMap->getGlobalElement(Bindices_i[k]);
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
    else {
      for(size_t k=0; k<Blen_i; ++k) {
        temp = bgids[Bindices_i[k]];
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
  }

  //Epetra_Util util;

  const bool C_filled = C.isFillComplete();

  //To form C = A*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(j,k) )
  //
  //This is the easiest case of all to code (easier than A*B, A^T*B, A^T*B^T).
  //But it requires the use of a 'sparsedot' function (we're simply forming
  //dot-products with row A_i and row B_j for all i and j).

  //loop over the rows of A.
  for(size_t i=0; i<Aview->numRows; ++i) {
    if (Aview->remote[i]) {
      continue;
    }

    ArrayView<const LocalOrdinal> Aindices_i = Aview->indices[i];
    ArrayView<const Scalar> Aval_i  = Aview->values[i];
    size_t A_len_i = Aview->numEntriesPerRow[i];

    if (A_len_i < 1) {
      continue;
    }

    for(size_t k=0; k<A_len_i; ++k) {
      Aind[k] = Aview->colMap->getGlobalElement(Aindices_i[k]);
      avals[k] = Aval_i[k];
    }


    typename ArrayView<GlobalOrdinal>::iterator end = Aind.begin() + A_len_i;
    //for(size_t l = 0; l<A_len_i; ++l, ++end);
    //util.Sort(true, A_len_i, Aind, 1, &avals, 0, NULL);
    sort2(Aind.begin(), end, avals.begin());

    //int mina = Aind[0];
    GlobalOrdinal mina = Aind[0];
    //int maxa = Aind[A_len_i-1];
    GlobalOrdinal maxa = Aind[A_len_i-1];

    if (mina > max_all_b || maxa < min_all_b) {
      continue;
    }

    GlobalOrdinal global_row = Aview->rowMap->getGlobalElement(i);

    //loop over the rows of B and form results C_ij = dot(A(i,:),B(j,:))
    for(size_t j=0; j<Bview->numRows; ++j) {
      if (b_firstcol[j] > maxa || b_lastcol[j] < mina) {
        continue;
      }

      ArrayView<const LocalOrdinal> Bindices_j = Bview->indices[j];

      size_t B_len_j = Bview->numEntriesPerRow[j];
      if (B_len_j < 1) {
        continue;
      }

      //int tmp, Blen = 0;
      GlobalOrdinal tmp;
      size_t Blen = 0;

      if (Bview->remote[j]) {
        for(size_t k=0; k<B_len_j; ++k) {
          tmp = Bview->importColMap->getGlobalElement(Bindices_j[k]);
          if (tmp < mina || tmp > maxa) {
            continue;
          }
          bvals[Blen] = Bview->values[j][k];
          Bind[Blen++] = tmp;
        }
      }
      else {
        for(size_t k=0; k<B_len_j; ++k) {
          tmp = bgids[Bindices_j[k]];
          if (tmp < mina || tmp > maxa) {
            continue;
          }
          bvals[Blen] = Bview->values[j][k];
          Bind[Blen++] = tmp;
        }
      }
      if (Blen < 1) {
        continue;
      }

      //util.Sort(true, Blen, Bind, 1, &bvals, 0, NULL);
      sort2(Bind.begin(), Bind.end(), bvals.begin());

      const Scalar C_ij = MMdetails::sparsedot<double,int>(avals(0,A_len_i), Aind(0,A_len_i), bvals(0,B_len_j), Bind(0,B_len_j));

      if (C_ij == ScalarTraits<Scalar>::zero()) {
        continue;
      }
      GlobalOrdinal global_col = Bview->rowMap->getGlobalElement(j);

      if (C_filled) {
        C.sumIntoGlobalValues(global_row, tuple(global_col), tuple(C_ij));
      }
      else {
        C.insertGlobalValues(global_row, tuple(global_col), tuple(C_ij));
      }

    }
  }

  return;
}


//kernel method for computing the local portion of C = A^T*B
template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::mult_Atrans_B(
  RCP<CrsMatrixStruct_t >& Aview, 
  RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>&  C)
{
  LocalOrdinal C_firstCol = Bview->colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview->colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = OrdinalTraits<LocalOrdinal>::invalid();

  if (Bview->importColMap != null) {
    C_firstCol_import = Bview->importColMap->getMinLocalIndex();
    C_lastCol_import = Bview->importColMap->getMaxLocalIndex();
  }

  size_t C_numCols = C_lastCol - C_firstCol + OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Array<Scalar> C_row_i = Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> C_colInds = Array<GlobalOrdinal>(C_numCols);

  size_t i, j, k;

  for(j=OrdinalTraits<size_t>::zero(); j<C_numCols; ++j) {
    C_row_i[j] = ScalarTraits<Scalar>::zero();
    C_colInds[j] = OrdinalTraits<GlobalOrdinal>::zero();
  }

  //To form C = A^T*B, compute a series of outer-product updates.
  //
  // for (ith column of A^T) { 
  //   C_i = outer product of A^T(:,i) and B(i,:)
  // Where C_i is the ith matrix update,
  //       A^T(:,i) is the ith column of A^T, and
  //       B(i,:) is the ith row of B.
  // }
  //

  //dumpCrsMatrixStruct(Aview->;
  //dumpCrsMatrixStruct(Bview);
  //int localProc = Bview.colMap->Comm().MyPID();
  int localProc = Bview->colMap->getComm()->getRank();

  ArrayView<const GlobalOrdinal> Arows = Aview->rowMap->getNodeElementList();

  bool C_filled = C.isFillComplete();

  //loop over the rows of A (which are the columns of A^T).
  for(i=OrdinalTraits<size_t>::zero(); i<Aview->numRows; ++i) {

    ArrayView<const LocalOrdinal> Aindices_i = Aview->indices[i];
    ArrayView<const Scalar> Aval_i  = Aview->values[i];

    //we'll need to get the row of B corresponding to Arows[i],
    //where Arows[i] is the GID of A's ith row.
    LocalOrdinal Bi = Bview->rowMap->getLocalElement(Arows[i]);
    TEST_FOR_EXCEPTION(Bi == OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
      "mult_Atrans_B ERROR, proc "<<localProc<<" needs row "
     <<Arows[i]<<" of matrix B, but doesn't have it.");

    ArrayView<const LocalOrdinal> Bcol_inds = Bview->indices[Bi];
    ArrayView<const Scalar> Bvals_i = Bview->values[Bi];

    //for each column-index Aj in the i-th row of A, we'll update
    //global-row GID(Aj) of the result matrix C. In that row of C,
    //we'll update column-indices given by the column-indices in the
    //ith row of B that we're now holding (Bcol_inds).

    //First create a list of GIDs for the column-indices
    //that we'll be updating.

    size_t Blen = Bview->numEntriesPerRow[Bi];
    if (Bview->remote[Bi]) {
      for(j=OrdinalTraits<size_t>::zero(); j<Blen; ++j) {
        C_colInds[j] = Bview->importColMap->getGlobalElement(Bcol_inds[j]);
      }
    }
    else {
      for(j=OrdinalTraits<size_t>::zero(); j<Blen; ++j) {
        C_colInds[j] = Bview->colMap->getGlobalElement(Bcol_inds[j]);
      }
    }

    //loop across the i-th row of A (column of A^T)
    for(j=OrdinalTraits<size_t>::zero(); j<Aview->numEntriesPerRow[i]; ++j) {

      LocalOrdinal Aj = Aindices_i[j];
      Scalar Aval = Aval_i[j];

      GlobalOrdinal global_row;
      if (Aview->remote[i]) {
        global_row = Aview->importColMap->getGlobalElement(Aj);
      }
      else {
        global_row = Aview->colMap->getGlobalElement(Aj);
      }

      if (!C.getRowMap()->isNodeGlobalElement(global_row)) {
        continue;
      }

      for(k=OrdinalTraits<size_t>::zero(); k<Blen; ++k) {
        C_row_i[k] = Aval*Bvals_i[k];
      }

      //
      //Now add this row-update to C.
      //

      C_filled ?
        C.sumIntoGlobalValues(global_row, C_colInds(), C_row_i() )
        :
        C.insertGlobalValues(global_row, C_colInds(), C_row_i());

    }
  }

  //delete [] C_row_i;
  //delete [] C_colInds;

  return(0);
}

//kernel method for computing the local portion of C = A^T*B^T
template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::mult_Atrans_Btrans(
  RCP<CrsMatrixStruct_t >& Aview, 
  RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{
  LocalOrdinal C_firstCol = Aview->rowMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Aview->rowMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = OrdinalTraits<LocalOrdinal>::invalid();

  if (Aview->importColMap != null) {
    C_firstCol_import = Aview->importColMap->getMinLocalIndex();
    C_lastCol_import = Aview->importColMap->getMaxLocalIndex();
  }

  size_t C_numCols = C_lastCol - C_firstCol + OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Array<Scalar> dwork = Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);

  Array<Scalar> C_col_j = dwork;
  Array<GlobalOrdinal> C_inds = iwork;

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);


  size_t i, j, k;

  for(j=OrdinalTraits<size_t>::zero(); j<C_numCols; ++j) {
    C_col_j[j] = ScalarTraits<Scalar>::zero();
    C_inds[j] = OrdinalTraits<GlobalOrdinal>::invalid();
  }

  ArrayView<const GlobalOrdinal> A_col_inds = Aview->colMap->getNodeElementList();
  ArrayView<const GlobalOrdinal> A_col_inds_import = Aview->importColMap == null ?
    Aview->importColMap->getNodeElementList() 
	:
	null;

  RCP<const Map_t > Crowmap = C.getRowMap();

  //To form C = A^T*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(k,i)*B(j,k) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc. In other words, we avoid
  //column-wise operations like the plague...

  ArrayView<const GlobalOrdinal> Brows = Bview->rowMap->getNodeElementList();

  //loop over the rows of B
  for(j=OrdinalTraits<size_t>::zero(); j<Bview->numRows; ++j) {
    ArrayView<const LocalOrdinal> Bindices_j = Bview->indices[j];
    ArrayView<const Scalar> Bvals_j = Bview->values[j];

    //GlobalOrdinal global_col = Brows[j];
    ArrayView<const GlobalOrdinal> global_col = Brows.view(j,1);

    //loop across columns in the j-th row of B and for each corresponding
    //row in A, loop across columns and accumulate product
    //A(k,i)*B(j,k) into our partial sum quantities in C_col_j. In other
    //words, as we stride across B(j,:), we use selected rows in A to
    //calculate updates for column j of the result matrix C.

    for(k=OrdinalTraits<size_t>::zero(); k<Bview->numEntriesPerRow[j]; ++k) {

      LocalOrdinal bk = Bindices_j[k];
      Scalar Bval = Bvals_j[k];

      GlobalOrdinal global_k;
      if (Bview->remote[j]) {
        global_k = Bview->importColMap->getGlobalElement(bk);
      }
      else {
        global_k = Bview->colMap->getGlobalElement(bk);
      }

      //get the corresponding row in A
      LocalOrdinal ak = Aview->rowMap->getLocalElement(global_k);
      if (ak == OrdinalTraits<LocalOrdinal>::invalid()) {
        continue;
      }

      ArrayView<const LocalOrdinal> Aindices_k = Aview->indices[ak];
      ArrayView<const Scalar> Avals_k = Aview->values[ak];

      size_t C_len = OrdinalTraits<size_t>::zero();

      if (Aview->remote[ak]) {
        for(i=OrdinalTraits<size_t>::zero(); i<Aview->numEntriesPerRow[ak]; ++i) {
          C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds_import[Aindices_k[i]];
        }
      }
      else {
        for(i=OrdinalTraits<size_t>::zero(); i<Aview->numEntriesPerRow[ak]; ++i) {
          C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds[Aindices_k[i]];
        }
      }

      //Now loop across the C_col_j values and put non-zeros into C.

      for(i=OrdinalTraits<size_t>::zero(); i < C_len ; ++i) {
        if (C_col_j[i] == ScalarTraits<Scalar>::zero()) continue;

        GlobalOrdinal global_row = C_inds[i];
        if (!Crowmap->isNodeGlobalElement(global_row)) {
          continue;
        }

/*  int err = C.SumIntoGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);

  if (err < 0) {
    return(err);
  }
  else {
          if (err > 0) {
      err = C.InsertGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);
      if (err < 0) {
              return(err);
            }
    }
  }*/

        try{
          //C.sumIntoGlobalValues(global_row, global_col, C_col_j[i]);
          C.sumIntoGlobalValues(global_row, global_col, C_col_j.view(i,1));
        }
        catch(std::runtime_error){
          //C.insertGlobalValues(global_row, global_col, C_col_j[i]);
          C.insertGlobalValues(global_row, global_col, C_col_j.view(i,1));
        }
      }
    }
  }

  //delete [] dwork;
  //delete [] iwork;

  return(0);
}

template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::import_and_extract_views(
  RCP<const CrsMatrix_t >& M,
  RCP<const Map_t >& targetMap,
  RCP<CrsMatrixStruct_t >& Mview)
{
  // CGB: this statement only returns 0. so, why bother with a return value? throw exception for error.

  // The goal of this method is to populate the 'Mview' struct with views of the
  // rows of M, including all rows that correspond to elements in 'targetMap'.
  // 
  // If targetMap includes local elements that correspond to remotely-owned rows
  // of M, then those remotely-owned rows will be imported into
  // 'Mview.importMatrix', and views of them will be included in 'Mview'.
  Mview->deleteContents();

  const RCP<const Map_t > Mrowmap = M->getRowMap();

  // CGB: will/should this change? then make it const, to ensure that it doesn't?
  //OLD// int numProcs = Mrowmap->getComm()->getSize();
  const int numProcs = Mrowmap->getComm()->getSize();

  ArrayView<const GlobalOrdinal> Mrows = targetMap->getNodeElementList();

  // CGB: Resize exists for a reason; use it
  //OLD//  Mview->numEntriesPerRow = Array<size_t>(Mview->numRows);
  //OLD//  Mview->indices = 
  //OLD//    Array<ArrayView<const LocalOrdinal> >(Mview->numRows);
  //OLD//  Mview->values = 
  //OLD//    Array<ArrayView<const Scalar> >(Mview->numRows);
  //OLD//  Mview->remote = Array<bool>(Mview->numRows);
  Mview->numRemote = 0;
  Mview->numRows = targetMap->getNodeNumElements();
  Mview->numEntriesPerRow.resize(Mview->numRows);
  Mview->indices.resize(         Mview->numRows);
  Mview->values.resize(          Mview->numRows);
  Mview->remote.resize(          Mview->numRows);
  Mview->origRowMap = M->getRowMap();
  Mview->rowMap = targetMap;
  Mview->colMap = M->getColMap();
  Mview->domainMap = M->getDomainMap();
  Mview->importColMap = null;

  // CGB: some documentation, e.g. of this loop, would be nice... Mike may not need it, but you will likely not be the one maintaining this code, so you should do a courtesy 
  // to those who will follow you

  // mark each row in targetMap as local or remote, and go ahead and get a view for the local rows

  // CGB: this is C++, not FORTRAN 66. If you don't need this outside of the loop, then don't pollute the local scope with it.
  //OLD// size_t i;
  for(size_t i=0; i < Mview->numRows; ++i) 
  {
    // CGB: is it gonna change? no? make it const
    //OLD// LocalOrdinal mlid = Mrowmap->getLocalElement(Mrows[i]);
    const LocalOrdinal mlid = Mrowmap->getLocalElement(Mrows[i]);
    // CGB: this is not correct; for non-locals, mlid will be invalid(), not negative one. /
    //OLD// if (mlid < OrdinalTraits<LocalOrdinal>::zero()) 
    // See http://trilinos.sandia.gov/packages/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a11c6a4585d616718eeb5c16d8a492511

    if (mlid == OrdinalTraits<LocalOrdinal>::invalid()) {
      Mview->remote[i] = true;
      ++Mview->numRemote;
    }
    else {
      Mview->remote[i] = false;
      M->getLocalRowView(mlid, Mview->indices[i], Mview->values[i]);
	    Mview->numEntriesPerRow[i] = Mview->indices[i].size();
    }
  }

  if (numProcs < 2) {
    TEST_FOR_EXCEPTION(Mview->numRemote > 0, std::runtime_error,
      "MatrixMatrix::import_and_extract_views ERROR, numProcs < 2 but attempting to import remote matrix rows." <<std::endl);
    // CGB: the exception handles the error; this return statement is erroneous
    // return(-1);

    //If only one processor we don't need to import any remote rows, so return.
    return(0);
  }

  //
  // Now we will import the needed remote rows of M, if the global maximum
  // value of numRemote is greater than 0.
  //

  global_size_t globalMaxNumRemote = 0;
  Teuchos::reduceAll(*(Mrowmap->getComm()) , Teuchos::REDUCE_MAX, Mview->numRemote, Teuchos::outArg(globalMaxNumRemote) );

  if (globalMaxNumRemote > 0) {

    // Create a map that describes the remote rows of M that we need.

    Array<GlobalOrdinal> MremoteRows(Mview->numRemote);

    // CGB: I know where this came from, but it doesn't need to be here. More danger from trying to change Epetra to Tpetra, one line at at time.
    //OLD// if(Mview->numRemote > 0) {
    //OLD//  Array<GlobalOrdinal>(Mview->numRemote);
    //OLD// }

    global_size_t offset = 0;
    for(size_t i=0; i < Mview->numRows; ++i) {
      if (Mview->remote[i]) {
        MremoteRows[offset++] = Mrows[i];
      }
    }

    RCP<Map_t > MremoteRowMap = rcp(new Map_t(OrdinalTraits<GlobalOrdinal>::invalid(), MremoteRows(), Mrowmap->getIndexBase(), Mrowmap->getComm(), Mrowmap->getNode()));

    // Create an importer with target-map MremoteRowMap and source-map Mrowmap.
    Import<LocalOrdinal, GlobalOrdinal, Node> importer(Mrowmap, MremoteRowMap);

    // Now create a new matrix into which we can import the remote rows of M that we need.
    Mview->importMatrix = rcp(new CrsMatrix<Scalar,LocalOrdinal, GlobalOrdinal, Node, SpMatOps>( MremoteRowMap, 1 ));
    Mview->importMatrix->doImport(*M, importer, INSERT);
    Mview->importMatrix->fillComplete(M->getDomainMap(), M->getRangeMap());

    // Save the column map of the imported matrix, so that we can convert indices back to global for arithmetic later
    Mview->importColMap = Mview->importMatrix->getColMap();

    // Finally, use the freshly imported data to fill in the gaps in our views of rows of M.
    for(size_t i=0; i < Mview->numRows; ++i) 
    {
      if (Mview->remote[i]) {
        const LocalOrdinal importLID = MremoteRowMap->getLocalElement(Mrows[i]);
        Mview->importMatrix->getLocalRowView(importLID,
                                             Mview->indices[i],
                                             Mview->values[i]);
        Mview->numEntriesPerRow[i] = Mview->indices[i].size();
      }
    }
  }
  return(0);
}

// CGB: check this...
template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
template<class Ordinal >
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::distribute_list(
  const RCP<const Teuchos::Comm<Ordinal> > comm,
  size_t lenSendList,
  const Array<GlobalOrdinal>& sendList,
  size_t& maxSendLen,
  Array<GlobalOrdinal>& recvList)
{
  maxSendLen = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, OrdinalTraits<Ordinal>::one(), &lenSendList, &maxSendLen);
  const int numProcs = comm->getSize();
  recvList.resize(numProcs*maxSendLen);
  Array<GlobalOrdinal> send(maxSendLen);
  std::copy(sendList.begin(), sendList.end(), send.begin());
  Teuchos::gatherAll(*comm, (Ordinal)maxSendLen, send.getRawPtr(), (Ordinal)(numProcs*maxSendLen), recvList.getRawPtr());
  return(0);
}


template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > 
MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::create_map_from_imported_rows(
  RCP<const Map_t > map,
  size_t totalNumSend,
  const ArrayView<const GlobalOrdinal> &sendRows,
  int numProcs,
  const ArrayView<const size_t> &numSendPerProc)
{
  // Perform sparse all-to-all communication to send the row-GIDs
  // in sendRows to appropriate processors according to offset
  // information in numSendPerProc.
  // Then create and return a map containing the rows that we
  // received on the local processor.

  RCP<Distributor> distributor = rcp(new Distributor(map->getComm()));
  Array<int> sendPIDs(totalNumSend);
  int offset = 0;
  for(int i=0; i<numProcs; ++i) {
    for(size_t j=0; j<numSendPerProc[i]; ++j) {
      sendPIDs[offset++] = i;
    }
  }

  size_t numRecv = 0;
  numRecv = distributor->createFromSends(sendPIDs());

  Array<GlobalOrdinal> recv_rows(numRecv);
  const size_t numpackets = 1;
  distributor->doPostsAndWaits(sendRows.getConst(), numpackets, recv_rows());

  //Now create a map with the rows we've received from other processors.
  RCP<Map_t > import_rows = rcp(new Map_t(OrdinalTraits<global_size_t>::invalid(), recv_rows(), map->getIndexBase(), map->getComm()));

  return( import_rows );
}


template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::form_map_union(RCP<const Map_t> map1, RCP<const Map_t> map2)
{
  // form the union of two maps
  if (map1 == null) {
    return map2;
  }
  if (map2 == null) {
    return map1;
  }

  const size_t map1_len = map1->getNodeNumElements();
  ArrayView<const GlobalOrdinal> map1_elements = map1->getNodeElementList();
  const size_t map2_len = map2->getNodeNumElements();
  ArrayView<const GlobalOrdinal> map2_elements = map2->getNodeElementList();

  Array<GlobalOrdinal> union_elements(map1_len + map2_len);
    
  //int map1_offset = 0, map2_offset = 0, union_offset = 0;
  size_t map1_offset = 0;
  size_t map2_offset = 0;
  size_t union_offset = 0;

  while(map1_offset < map1_len && map2_offset < map2_len) {
    GlobalOrdinal map1_elem = map1_elements[map1_offset];
    GlobalOrdinal map2_elem = map2_elements[map2_offset];
    if (map1_elem < map2_elem) {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
    }
    else if (map1_elem > map2_elem) {
      union_elements[union_offset++] = map2_elem;
      ++map2_offset;
    }
    else {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
      ++map2_offset;
    }
  }
  for(size_t i=map1_offset; i<map1_len; ++i) {
    union_elements[union_offset++] = map1_elements[i];
  }
  for(size_t i=map2_offset; i<map2_len; ++i) {
    union_elements[union_offset++] = map2_elements[i];
  }

  // note, union_elements potentially contains duplicate elements. however, this is 
  // allowed, and accounted for, in the Map constructor below.

  // CGB: I think this should be union_offset, not +1
  // ERROR FIXED
  //OLD// union_elements.resize(union_offset+1);
  union_elements.resize(union_offset);
  return rcp(new Map_t(OrdinalTraits<global_size_t>::invalid(), union_elements(), map1->getIndexBase(), map1->getComm(), map1->getNode()));
}


template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps >
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::find_rows_containing_cols(
  RCP<const CrsMatrix_t > M,
  RCP<const Map_t > colmap)
{
  //The goal of this function is to find all rows in the matrix M that contain
  //column-indices which are in 'colmap'. A map containing those rows is
  //returned.

  const int numProcs = colmap->getComm()->getSize();
  const int localProc = colmap->getComm()->getRank();

  if (numProcs < 2) {
    return M->getRowMap();
  }

  const size_t MnumRows = M->getNodeNumRows();
  const size_t numCols = colmap->getNodeNumElements();

  Array<GlobalOrdinal> cols(numCols + 1);
  cols[0] = numCols;
  // get ids from column map and sort them
  cols(1,numCols).assign(colmap->getNodeElementList());
  // must sort before distribute_list, because below we assume that all ids are sorted for each proc
  std::sort(cols.begin()+1, cols.end());

  size_t max_num_cols;
  Array<GlobalOrdinal> all_proc_cols;
  distribute_list(colmap->getComm(), numCols+1, cols, max_num_cols, all_proc_cols);

  RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > Mgraph = M->getCrsGraph();
  RCP<const Map_t > Mrowmap = M->getRowMap();
  RCP<const Map_t > Mcolmap = M->getColMap();

  Array<size_t> procNumCols(numProcs);
  Array<size_t> procNumRows(numProcs);
  Array<GlobalOrdinal> procRows_1D(numProcs*MnumRows);
  Array<typename Array<GlobalOrdinal>::iterator> procCols(numProcs);
  Array<typename Array<GlobalOrdinal>::iterator> procRows(numProcs);
  size_t offset = 0;
  for(int p=0; p<numProcs; ++p) {
    // procCols points into all_proc_cols, and procNumCols[p] stores the number of cols for proc p
    procNumCols[p] = all_proc_cols[offset];
    procCols[p] = all_proc_cols.begin() + offset + 1;
    offset += max_num_cols;
    // procRows points into procRows_1D, and procNumRows[p] stores ??? FINISH
    procNumRows[p] = 0;
    procRows[p] = procRows_1D.begin() + p*MnumRows;
  }


  for(LocalOrdinal localRow =  Mrowmap->getMinLocalIndex(); 
                   localRow <= Mrowmap->getMaxLocalIndex(); 
                 ++localRow) 
  {
    const GlobalOrdinal globalRow = Mrowmap->getGlobalElement(localRow);
    ArrayView<const LocalOrdinal> Mindices;
    Mgraph->getLocalRowView(localRow, Mindices);
    for (size_t j=0; j<(size_t)Mindices.size(); ++j) {
      const GlobalOrdinal colGID = Mcolmap->getGlobalElement(Mindices[j]);
      for(int p=0; p<numProcs; ++p) 
      {
        if (p==localProc) continue;
        // according to the sort above, before distribute_list, these are sorted
        bool result = std::binary_search(procCols[p], procCols[p] + procNumCols[p], colGID);
        if (result) {
          size_t numRowsP = procNumRows[p];
          typename ArrayView<GlobalOrdinal>::iterator prows = procRows[p];
          if (numRowsP < 1 || prows[numRowsP-1] < globalRow) {
            prows[numRowsP] = globalRow;
            procNumRows[p]++;
          }
        }
      }
    }
  }

  // Now make the contents of procRows occupy a contiguous section
  // of procRows_1D; pack the last numProcs-1 sections.
  offset = procNumRows[0];
  for(int i=1; i<numProcs; ++i) {
    // we are packing, which must be done sequentially
    for(size_t j=0; j<procNumRows[i]; ++j) {
      procRows_1D[offset++] = procRows[i][j];
    }
  }
  const size_t totalNumSend = offset;

  // Next we will do a sparse all-to-all communication to send the lists of rows
  // to the appropriate processors, and create a map with the rows we've received
  // from other processors.
  RCP<const Map_t > recvd_rows = 
    create_map_from_imported_rows(
      Mrowmap,
      totalNumSend,
      procRows_1D(),
      numProcs,
      procNumRows());
  RCP<const Map_t > result_map = form_map_union(M->getRowMap(), recvd_rows);
  return(result_map);
}


// CGB: check this...
template <class Scalar, 
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class SpMatOps >
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::Multiply(
  RCP<const CrsMatrix_t > A,
  bool transposeA,
  RCP<const CrsMatrix_t > B,
  bool transposeB,
  RCP<CrsMatrix_t > C,
  bool call_FillComplete_on_result)
{
  //
  //This method forms the matrix-matrix product C = op(A) * op(B), where
  //op(A) == A   if transposeA is false,
  //op(A) == A^T if transposeA is true,
  //and similarly for op(B).
  //

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix A is not fill complete.");
  TEST_FOR_EXCEPTION(!B->isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix B is not fill complete.");

  //We're going to refer to the different combinations of op(A) and op(B)
  //as scenario 1 through 4.

  int scenario = 1;//A*B
  if (transposeB && !transposeA) scenario = 2;//A*B^T
  if (transposeA && !transposeB) scenario = 3;//A^T*B
  if (transposeA && transposeB)  scenario = 4;//A^T*B^T

  //now check size compatibility
  global_size_t Aouter = transposeA ? A->getGlobalNumCols() : A->getGlobalNumRows();
  global_size_t Bouter = transposeB ? B->getGlobalNumRows() : B->getGlobalNumCols();
  global_size_t Ainner = transposeA ? A->getGlobalNumRows() : A->getGlobalNumCols();
  global_size_t Binner = transposeB ? B->getGlobalNumCols() : B->getGlobalNumRows();
  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    <<Aouter<<"x"<<Ainner << ", op(B) is "<<Binner<<"x"<<Bouter<<std::endl);

  //The result matrix C must at least have a row-map that reflects the
  //correct row-size. Don't check the number of columns because rectangular
  //matrices which were constructed with only one map can still end up
  //having the correct capacity and dimensions when filled.
  TEST_FOR_EXCEPTION(Aouter > C->getGlobalNumRows(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has "<<C->getGlobalNumRows()
     << " rows, should have at least "<<Aouter << std::endl);

  //It doesn't matter whether C is already Filled or not. If it is already
  //Filled, it must have space allocated for the positions that will be
  //referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  //we'll error out later when trying to store result values.
  
  // CGB: However, matrix must be in active-fill
  TEST_FOR_EXCEPT( C->isFillActive() == false );

  //We're going to need to import remotely-owned sections of A and/or B
  //if more than 1 processor is performing this run, depending on the scenario.
  int numProcs = A->getComm()->getSize();

  //If we are to use the transpose of A and/or B, we'll need to be able to 
  //access, on the local processor, all rows that contain column-indices in
  //the domain-map.
//  const Epetra_Map* domainMap_A = &(A.DomainMap());
//  const Epetra_Map* domainMap_B = &(B.DomainMap());

  //const Epetra_Map* rowmap_A = &(A.RowMap());
  //const Epetra_Map* rowmap_B = &(B.RowMap());

  //Declare some 'work-space' maps which may be created depending on
  //the scenario, and which will be deleted before exiting this function.


  //Declare a couple of structs that will be used to hold views of the data
  //of A and B, to be used for fast access during the matrix-multiplication.
  RCP<CrsMatrixStruct_t > Aview = rcp(new CrsMatrixStruct_t);
  RCP<CrsMatrixStruct_t > Bview = rcp(new CrsMatrixStruct_t);

  //const Epetra_Map* targetMap_A = rowmap_A;
  //const Epetra_Map* targetMap_B = rowmap_B;

  RCP<const Map_t > targetMap_A = A->getRowMap();
  RCP<const Map_t > targetMap_B = B->getRowMap();

  if (numProcs > 1) {
    //If op(A) = A^T, find all rows of A that contain column-indices in the
    //local portion of the domain-map. (We'll import any remote rows
    //that fit this criteria onto the local processor.)
    if (transposeA) {
      targetMap_A = find_rows_containing_cols(A, A->getDomainMap());
    }
  }
  //Now import any needed remote rows and populate the Aview struct.
  //EPETRA_CHK_ERR( import_and_extract_views(A, *targetMap_A, Aview) );
  import_and_extract_views(A, targetMap_A, Aview);

  //We will also need local access to all rows of B that correspond to the
  //column-map of op(A).
  if (numProcs > 1) {
    //const Epetra_Map* colmap_op_A = NULL;
    RCP<const Map_t > colmap_op_A = null;
    if (transposeA) {
      colmap_op_A = targetMap_A;
    }
    else {
      colmap_op_A = A->getColMap(); 
    }

    targetMap_B = colmap_op_A;

    //If op(B) = B^T, find all rows of B that contain column-indices in the
    //local-portion of the domain-map, or in the column-map of op(A).
    //We'll import any remote rows that fit this criteria onto the
    //local processor.
    if (transposeB) {
      RCP<const Map_t > mapunion1 = form_map_union(colmap_op_A, B->getDomainMap());
      if (MMdebug::debug_level != Teuchos::VERB_NONE) {
        *MMdebug::debug_stream << "mapunion1" << std::endl;
        mapunion1->describe(*MMdebug::debug_stream, MMdebug::debug_level);
      }
      targetMap_B = find_rows_containing_cols(B, mapunion1);
    }
  }
  if (MMdebug::debug_level != Teuchos::VERB_NONE) {
    *MMdebug::debug_stream << "targetMap_B" << std::endl;
    targetMap_B->describe(*MMdebug::debug_stream, MMdebug::debug_level);
  }

  //Now import any needed remote rows and populate the Bview struct.
  import_and_extract_views(B, targetMap_B, Bview);

  if (MMdebug::debug_level != Teuchos::VERB_NONE) {
    *MMdebug::debug_stream << "C->getRowMap()" << std::endl;
    C->getRowMap()->describe(*MMdebug::debug_stream, MMdebug::debug_level);
  }

  //If the result matrix C is not already FillComplete'd, we will do a
  //preprocessing step to create the nonzero structure, then call FillComplete,
  if (!C->isFillComplete()) {
    CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsgraphbuilder(C->getRowMap());

    //pass the graph-builder object to the multiplication kernel to fill in all
    //the nonzero positions that will be used in the result matrix.
    switch(scenario) {
    case 1:    mult_A_B(Aview, Bview, crsgraphbuilder);
      break;
    case 2:    mult_A_Btrans(Aview, Bview, crsgraphbuilder);
      break;
    case 3:    mult_Atrans_B(Aview, Bview, crsgraphbuilder);
      break;
    case 4:    mult_Atrans_Btrans(Aview, Bview, crsgraphbuilder);
      break;
    }

    //now insert all of the nonzero positions into the result matrix.
    insert_matrix_locations(crsgraphbuilder, C);

  /*  if (call_FillComplete_on_result) {
      RCP<const Map_t > domainmap = transposeB ? B->getRangeMap() : B->getDomainMap();

      RCP<const Map_t > rangemap = transposeA ? A->getDomainMap() : A->getRangeMap();

      C->fillComplete(domainmap, rangemap, DoNotOptimizeStorage);
      call_FillComplete_on_result = false;
    }*/
  }
  //RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //C->describe(*out, Teuchos::VERB_EXTREME);

  //Pre-zero the result matrix:
  //C->setAllToScalar(ScalarTraits<Scalar>::zero());

  //Now call the appropriate method to perform the actual multiplication.

  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  switch(scenario) {
  case 1:    mult_A_B(Aview, Bview, crsmat);
    break;
  case 2:    mult_A_Btrans(Aview, Bview, crsmat);
    break;
  case 3:    mult_Atrans_B(Aview, Bview, crsmat);
    break;
  case 4:    mult_Atrans_Btrans(Aview, Bview, crsmat);
    break;
  }

  if (call_FillComplete_on_result) {
    //We'll call FillComplete on the C matrix before we exit, and give
    //it a domain-map and a range-map.
    //The domain-map will be the domain-map of B, unless
    //op(B)==transpose(B), in which case the range-map of B will be used.
    //The range-map will be the range-map of A, unless
    //op(A)==transpose(A), in which case the domain-map of A will be used.
    if (!C->isFillComplete()) {
      RCP<const Map_t > domainmap = transposeB ? B->getRangeMap() : B->getDomainMap();

      RCP<const Map_t > rangemap = transposeA ? A->getDomainMap() : A->getRangeMap();
      //C->fillComplete(transposeB ? B->getRangeMap() : B->getDomainMap(), transposeA ? A->getDomainMap() : B->getRangeMap());
      C->fillComplete(domainmap, rangemap);
    }
  }

  return(0);
}

// CGB: check this...
template <class Scalar, 
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps >
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::Add(
  RCP<CrsMatrix_t > A,
  bool transposeA,
  Scalar scalarA,
  RCP<CrsMatrix_t > B,
  Scalar scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  //A should already be Filled. It doesn't matter whether B is
  //already Filled, but if it is, then its graph must already contain
  //all nonzero locations that will be referenced in forming the
  //sum.

  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix A.isFillComplete() is false, it is required to be true. (Result matrix B is not required to be isFillComplete()).");

  //explicit tranpose A formed as necessary
  RCP<CrsMatrix_t > Aprime = null;
  if( transposeA )
  {
	RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer(A);
    theTransposer.createTranspose(DoOptimizeStorage, Aprime);
  }
  else{
    Aprime = Teuchos::rcp_const_cast<CrsMatrix_t >(A);
  }

  /*size_t MaxNumEntries = std::max(
    A->getNodeMaxNumRowEntries(), 
    B->getNodeMaxNumRowEntries() );*/
  //int A_NumEntries, B_NumEntries;
  //size_t A_NumEntries, B_NumEntries;
  //int * A_Indices = new int[MaxNumEntries];
  size_t A_NumEntries;
  Array<GlobalOrdinal> A_Indices(Aprime->getGlobalMaxNumRowEntries());
  //double * A_Values = new double[MaxNumEntries];
  Array<Scalar> A_Values(Aprime->getGlobalMaxNumRowEntries());
  //int* B_Indices;
  //double* B_Values;
  ArrayView<const GlobalOrdinal> B_Indices;
  ArrayView<const Scalar> B_Values;

  size_t NumMyRows = B->getNodeNumRows();
  GlobalOrdinal Row;
  int ierr = 0;

  if( scalarB != ScalarTraits<Scalar>::zero() &&
    scalarB != ScalarTraits<Scalar>::one())
  {
    B->scale(scalarB);
  }
  
  if( scalarA )
  {
    //Loop over B's rows and sum into
    for( 
      size_t i = OrdinalTraits<size_t>::zero(); 
      i < NumMyRows;
      ++i )
    {
	    Row = B->getRowMap()->getGlobalElement(i);
      A_NumEntries = Aprime->getNumEntriesInGlobalRow(Row);
	    Aprime->getGlobalRowCopy(Row, A_Indices(), A_Values(), A_NumEntries);
      if (scalarA != ScalarTraits<Scalar>::one()) {
        for( 
          size_t j = OrdinalTraits<size_t>::zero(); 
          j < A_NumEntries; 
          ++j ) 
        {
          A_Values[j] *= scalarA;
        }
      
        if( B->isFillComplete() ) {//Sum In Values
          B->sumIntoGlobalValues( Row, A_Indices, A_Values );
        }
        else {
          B->insertGlobalValues( Row, A_Indices, A_Values);
        }
      }
    }
  }
  else {
      B->scale(scalarB);
  }
  return(ierr);
}

// CGB: check this...
template <class Scalar, 
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps>
int MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::Add(
  RCP<const CrsMatrix_t > A,
  bool transposeA,
  Scalar scalarA,
  RCP<const CrsMatrix_t > B,
  bool transposeB,
  Scalar scalarB,
  RCP<CrsMatrix_t > C)
{
  //
  //This method forms the matrix-matrix sum C = scalarA * op(A) + scalarB * op(B), where

  //A and B should already be Filled. C should be an empty pointer.


  TEST_FOR_EXCEPTION(!A.isFillComplete() || !B.isFillComplete(), std::runtime_error,
    "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() or B.Filled() is false,"
    "they are required to be true. (Result matrix C should be an empty pointer)" << std::endl);


  RCP<CrsMatrix_t > Aprime = null;
  RCP<CrsMatrix_t > Bprime = null;


  //explicit tranpose A formed as necessary
  if( transposeA ) {
	RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer(A);
    theTransposer.createTranspose(DoOptimizeStorage, Aprime);
  }
  else{
    Aprime = Teuchos::rcp_const_cast<CrsMatrix_t >(A);
  }

  //explicit tranpose B formed as necessary
  if( transposeB ) {
	RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer(B);
    theTransposer.createTranspose(DoOptimizeStorage, Bprime);
  }
  else{
    Bprime = Teuchos::rcp_const_cast<CrsMatrix_t >(B);
  }

  // allocate or zero the new matrix
  if(!C.is_null())
     C->setAllToScalar(ScalarTraits<Scalar>::zero());
  else
     C = rcp(new CrsMatrix_t(Aprime->getRowMap(), null));

  // build arrays  for easy resuse
  int ierr = 0;
  //Epetra_CrsMatrix * Mat[] = { Aprime,Bprime};
  Array<CrsMatrix_t > Mat = 
    Teuchos::tuple<CrsMatrix_t >(Aprime, Bprime);
  //double scalar[] = { scalarA, scalarB};
  Array<Scalar> scalar = Teuchos::tuple<Scalar>(scalarA, scalarB);

  // do a loop over each matrix to add: A reordering might be more efficient
  for(int k=0;k<2;k++) {
    size_t MaxNumEntries = Mat[k]->getNodeMaxNumRowEntries();
    size_t NumEntries;
    Array<LocalOrdinal> Indices;
    Array<Scalar> Values;
   
     size_t NumMyRows = Mat[k]->getNodeNumRows();
     GlobalOrdinal Row;
     int ierr = 0;
   
     //Loop over rows and sum into C
     for( size_t i = OrdinalTraits<size_t>::zero(); i < NumMyRows; ++i ) {
        Row = Mat[k]->getRowMap()->getGlobalElement(i);
		Mat[k]->extractGlobalRowCopy(Row, Indices, Values);
		NumEntries = Indices.size();
   
        if( scalar[k] != ScalarTraits<Scalar>::one() )
           for( size_t j = OrdinalTraits<size_t>::zero(); j < NumEntries; ++j ) Values[j] *= scalar[k];
   
        if(C->isFillComplete()) { // Sum in values
           C->sumIntoGlobalValues( Row, Indices, Values);
        } else { // just add it to the unfilled CRS Matrix
           C->insertGlobalValues( Row, Indices, Values);
        }
     }
  }
  return(ierr);
}


}

template<class Scalar, class LocalOrdinal>
Scalar Tpetra::MMdetails::sparsedot(
          const ArrayView<const Scalar>& u, const ArrayView<const LocalOrdinal>& u_ind,
          const ArrayView<const Scalar>& v, const ArrayView<const LocalOrdinal>& v_ind)
{
  const size_t usize = (size_t)u.size();
  const size_t vsize = (size_t)v.size();
  Scalar result = ScalarTraits<Scalar>::zero();
  size_t v_idx = 0;
  size_t u_idx = 0;
  while(v_idx < vsize && u_idx < usize) {
    LocalOrdinal ui = u_ind[u_idx];
    LocalOrdinal vi = v_ind[v_idx];
    if (ui < vi) {
      ++u_idx;
    }
    else if (ui > vi) {
      ++v_idx;
    }
    else {
      result += u[u_idx++]*v[v_idx++];
    }
  }
  return(result);
}


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MATRIXMATRIX_INSTANT(SCALAR,LO,GO,NODE,SPMATOPS) \
  \
  template<> \
  int MatrixMatrix::Multiply( \
    RCP<const CrsMatrix< SCALAR , LO , GO , NODE , SPMATOPS > >& A, \
    bool transposeA, \
    RCP<const CrsMatrix< SCALAR , LO , GO , NODE , SPMATOPS > >& B, \
    bool transposeB, \
    RCP<CrsMatrix< SCALAR , LO , GO , NODE , SPMATOPS > >& C, \
    bool call_FillComplete_on_result) \
  \
  template <> \
  int MatrixMatrix::Add( \
  RCP<const CrsMatrix< SCALAR, LO , GO , NODE ,  SPMATOPS > > A, \
  bool transposeA, \
  double scalarA, \
  RCP<CrsMatrix< SCALAR, LO , GO , NODE ,  SPMATOPS > > B, \
  double scalarB ) \
\
  template <> \
  int MatrixMatrix::Add( \
    RCP<const CrsMatrix< SCALAR, LO , GO , NODE ,  SPMATOPS > > A, \
    bool transposeA, \
    Scalar scalarA, \
    RCP<const CrsMatrix< SCALAR, LO , GO , NODE ,  SPMATOPS > > B, \
    bool transposeB, \
    Scalar scalarB, \
    RCP<CrsMatrix< SCALAR, LO , GO , NODE ,  SPMATOPS > > C) \
\
  template<> \
  RCP<const Map< LO , GO , NODE > > \
  MatrixMatrix::find_rows_containing_cols( \
    RCP<const CrsMatrix< SCALAR , LO , GO , NODE ,  SPMATOPS > > M, \
    RCP<const Map< LO , GO , NODE > > colmap)

#endif // TPETRA_MATRIXMATRIX_DEF_HPP
