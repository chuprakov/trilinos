// ///////////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_rank_2_chol_update.hpp
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

#ifndef RANK_2_CHOL_UPDATE_H
#define RANK_2_CHOL_UPDATE_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Update the cholesky factor of symmetric positive definite matrix for a rank-2 update.
 *
 * This function updates and an upper or lower triangular cholesky
 * factor in O(n^2) flops.  The basic idea is that the original symmetric
 * positive definite matrix
 \verbatim

       B = op(R')*op(R)
 \endverbatim
 * is subjected to a rank-2 update.  This rank-2 update to B must be expressed
 * in the form
 \verbatim

       B_new = J'*J

   where:
       J = (op(R) + a*u*v')
 \endverbatim
 * The basic approach is to introduce an othogonal matrix Q such that
 \verbatim

       Q*(op(R) + a*u*v') -> op(R_new)
 \endverbatim
 * is triangular.  In other words (with I = Q'*Q)
 \verbatim

       B_new = (op(R') + a*v*u')*Q'*Q*(op(R) + a*u*v') = op(R_new')*op(R_new)
 \endverbatim
 * This algorithm is based on Algorithm A3.4.1a in Dennis and Schnabel.
 * However, it has been extended in the sense that it handles upper and
 * lower along with transposed and nontransposed factors.  It does not
 * touch any of the elements of R outside of the designated triangular
 * region and returns the rotations used with no extra storage space
 * needed.  This is a much nicer function than any other one that I have
 * seen in the public domain.
 *  
 * It is guarrenteed that after this function finishes that the diagonal
 * of R_new will be positive.
 *
 * @param  a  [in] Scalar for update (see above).
 * @param  u  [in/out] (size n ) On input, contaitns the update vector u
 *            (see above).  On output, contains the rotations needed
 *            for the transformation Q1*u -> ||u||2 * e(last_i).
 *            (ToDo: Document this!).
 * @param  v  [in] (size n )Update vector (see above).
 * @param  w  [out] (size n-1) On output, contains stored rotations from
 *            BLAS_Cpp::rotg(...) used to transform the upper or lower
 *            Hessenberg itermediate factor R2 back to upper or lower
 *            triangular form.  (ToDo: Document this!).  If n == 1 then
 *            w can be NULL.
 * @param  R  [in/out] (size n x n) Upper or lower triangular factor.
 *            On input the diagonal must contain all positive elements on
 *            the diagonal.  On output contains the updated factor with
 *            all positive diagonal elements.
 * @param  R_trans
 *            [in] Determines if op(R) = R (no_trans) or op(R) = R' (trans)
 */
void rank_2_chol_update(
	const value_type     a
	,DVectorSlice         *u
	,const DVectorSlice   &v
	,DVectorSlice         *w
	,DMatrixSliceTriEle         *R
	,BLAS_Cpp::Transp    R_trans
	);

}  // end namespace AbstractLinAlgPack

#endif // RANK_2_CHOL_UPDATE_H
