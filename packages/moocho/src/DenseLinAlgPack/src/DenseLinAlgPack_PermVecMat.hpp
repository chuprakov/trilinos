// ///////////////////////////////////////////////////////////////////////////////////////
// PermVecMat.h
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
  * old vector.  After the permutation is perform the postcondition:\\
  * #     perm_vs(i) == vs(perm(i))#, for i = 1,2,...,#vs.size()#\\
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
  * Preconditions: \begin{itemize}
  * \item #perm.size() > 0# (throw std::length_error)
  * \end{itemize}
  *
  * Postconditions: \begin{itemize}
  * \item #perm(i) = i#, for i = 1, 2, ... m #perm.size()#
  * \end{itemize}
  */
void identity_perm(IVector* perm);

///
/** Find the inverse permutation (inv_perm) given a permutation vector (perm).
  *
  * Postconditions: \begin{itemize}
  * \item #inv_perm.size() == perm.size()#
  * \end{itemize}
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

#endif