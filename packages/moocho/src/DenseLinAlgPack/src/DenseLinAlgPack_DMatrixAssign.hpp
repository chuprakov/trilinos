// //////////////////////////////////////////////////////////////////////////////////
// GenMatrixAssign.h
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
// These are assignment functions that used to be in GenMatrixOp but because of
// some circular dependency problems their declarations where moved into a seprate
// file.

#ifndef GEN_MATRIX_ASSIGN_H
#define GEN_MATRIX_ASSIGN_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

// ///////////////////////////////////////////////////////////////////////////////////
/** @name {\bf GenMatrix/GenMatrixSlice Assignment Fucntions}.
  *
  * These are functions for assigning one matrix to another.  The rhs matrices
  * can be transposed (op(rhs) = rhs) or non-transposed (op(rhs) = trans(rhs).
  * The assignment operators for
  * GenMatrix and GenMatrixSlice use these functions to implement their
  * assignment operator functions.  The assignmet functions for triangular
  * matrices are ment to be called from the assignment operator functions
  * for wrapper classes for triangular (upper and lower) and sysmetric 
  * (upper and lower storage).
  *
  * These assignment functions check for aliasing and assignment to self
  * and create any needed temporaries if needed.
  */
//@{

/// gm_lhs = alpha (elementwise)
void assign(GenMatrix* gm_lhs, value_type alpha);

/// gm_lhs = op(gms_rhs)
void assign(GenMatrix* gm_lhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// gms_lhs = alpha (elementwise)
void assign(GenMatrixSlice* gms_lhs, value_type alpha);

/// gms_lhs = op(gms_rhs)
void assign(GenMatrixSlice* gms_lhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// tri_ele_gms_lhs = alpha (elementwise)
void assign(tri_ele_gms* tri_gms_lhs, value_type alpha);

/// tri_ele_gms_lhs = tri_ele_gms_rhs
void assign(tri_ele_gms* tri_ele_gms_lhs, const tri_ele_gms& tri_ele_gms_rhs);

//		end Assignment Fucntions
//@}

}	// end namespace LinAlgPack

#endif	// GEN_MATRIX_ASSIGN_H
