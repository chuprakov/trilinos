// ///////////////////////////////////////////////////////////////////////////////////////
// PermVecMat.h
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//
// Utilities for pivoting VectorSlices and GenMatrixSlices.

#ifndef PIVOTVECMAT_H
#define PIVOTVECMAT_H

#include <stdexcept>

#include "LinAlgPackTypes.h"

/** @name {\bf Vector / Matrix Permutations}.
  *
  * These are functions for pivoting the elements of a vector and the 
  * rows and/or columns of a rectandular matrix.
  *
  * For Vector and VectorSlice pivot funcitions the #perm# argument
  * gives the mapping from the new sequence to the old sequence.
  * The statement #i_old = perm(i_new)# gives the indice in the
  * old vector.  After the permutation is performed the postcondition:
  \verbatim

    perm_vs(i) == vs(perm(i)), for i = 1,2,...,vs.size()

  \endverbatim
  * is true.
  *
  * For the GenMatrix permutation functions the #row_perm# argument
  * gives the row permutations and the #col_perm# argument gives the column
  * permutations respectively.
  *
  * In each of these functions the permuted object can be the same
  * as the unpermuted object.
  */

//@{

//@Include: IVector.h

// /////////////////////////////////////////////////////////////////////////////////////////
// Public Permutation functions

namespace LinAlgPack {

///
/** Initialize a permutation to the identiy permutation.
  *
  * Preconditions: <ul>
  * <li> #perm.size() > 0# (throw std::length_error)
  * </ul>
  *
  * Postconditions: <ul>
  * <li> #perm(i) = i#, for i = 1, 2, ... m #perm.size()#
  * </ul>
  */
void identity_perm(IVector* perm);

///
/** Find the inverse permutation (inv_perm) given a permutation vector (perm).
  *
  * Postconditions: <ul>
  * <li> #inv_perm.size() == perm.size()#
  * </ul>
  */
void inv_perm(const IVector& perm, IVector* inv_perm);

/// Permute a VectorSlice in place
void perm_ele(const IVector& perm, VectorSlice* vs);

/// Perform y = P*x
void perm_ele(const VectorSlice& x, const IVector& perm, VectorSlice* y);

/// Perform x = P'*y
void inv_perm_ele(const VectorSlice& y, const IVector& perm, VectorSlice* x);

/// Permute a GenMatrixSlices rows
void perm_rows(const IVector& row_perm, GenMatrixSlice* gms);

/// Permute a GenMatrixSlices columns
void perm_cols(const IVector& col_perm, GenMatrixSlice* gms);

/// Permute a GenMatrixSlices rows and columns
void perm_rows_cols(const IVector& row_perm, const IVector& col_perm, GenMatrixSlice* gms);

} // end namespace LinAlgPack

//@}

#endif // PIVOTVECMAT_H
