// /////////////////////////////////////////////////////////////////
// DenseLinAlgPack_InvCholUpdate.hpp
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

#ifndef INV_CHOL_UPDATE_H
#define INV_CHOL_UPDATE_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

///
/* * Perform an update of the Cholesky factor of a symetric matrix (SM = L * L^T)
  * and while maintaining the triangular structure.
  *
  * This function permforms an update of L^T or L^-1 depending on whether
  * we are updating the cholesky factor or its inverse.  This function
  * implements algorithm A3.4.1 in Dennis and Schabel.  The update is:
  * (J_new^T = L_old^T + u * v^T) where J_new is rotated back to triangular form.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() == u.size() == v.size() (throw std::length_error)
  * </ul>
  */
void update_chol_factor(DMatrixSlice* UpTriM, DVectorSlice* u
	, const DVectorSlice& v);

///
/* * Perform a jacobi rotation or a matrix about row i.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() (throw std::length_error)
  * </ul>
  */
void jacobi_rotate(DMatrixSlice* UpTriM, size_type row_i, value_type alpha
	, value_type beta); 

}	// end namespace DenseLinAlgPack

#endif	// INV_CHOL_UPDATE_H
