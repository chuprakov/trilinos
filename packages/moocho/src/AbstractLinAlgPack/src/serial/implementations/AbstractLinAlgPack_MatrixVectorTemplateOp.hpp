// ///////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixVectorTemplateOp.hpp
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

#ifndef MATRIX_VECTOR_TEMPLATE_OP_H
#define MATRIX_VECTOR_TEMPLATE_OP_H

#include <stdexcept>
#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** @name {\bf Templated Matrix-DVector Operations}.
  *
  * These are the declarations (AbstractLinAlgPack_MatrixVectorTemplateOp.hpp) for template functions for
  * performing selected matrix-vector and matrix-matrix operations.  The templated
  * matrix type must have the following interface:
  * \begin{itemize}
  * \item #T_Matrix::size_type# - Member type for the type returned by #rows# and #cols#
  * \item #T_Matrix::row_type# - Member type for the type returned by #row(i)#
  * \item #T_Matrix::col_type# - Member type for the type returned by #col(j)#
  * \item #T_Matrix::size_type T_Matrix::rows() const#
  * \item #T_Matrix::size_type T_Matrix::cols() const#
  * \item #const T_Matrix::row_type& T_Matrix::row(size_type i) const# \
  *		- Returns a reference to the ith row (1,2,...,#T_Matrix::rows()#).
  * \item #const T_Matrix::col_type& T_Matrix::col(size_type j) const# \ 
  *		- Returns a reference to the ith row (1,2,...,#T_Rect_Matrix::rows()#).
  * \end{itemize}
  *
  * In addition, the following operations (UML notation) must be accessible for the types
  * #T_Matrix::row_type# and #T_Matrix::col_type# when these template functions
  * are instantiated:
  *
  * assign(v_lhs:DVector&,v_rhs:T_Matrix::const row_type&)
  * assign(vs_lhs:DVectorSlice&,v_rhs:const T_Matrix::row_type&)
  * dot(vs_rhs1:const DVectorSlice&, vs_rhs2: const T_Matrix::row_type&):value_type 
  * dot(vs_rhs1:const DVectorSlice&, vs_rhs2: const T_Matrix::col_type&):value_type 
  */
//@{

// ////////////////////////
// assign to a dense matrix

/// gm_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void assign(DMatrix& gm_lhs, const T_Matrix& gm_rhs, BLAS_Cpp::Transp trans_rhs);

/// gms_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void assign(DMatrixSlice& gms_lhs, const T_Matrix& gm_rhs, BLAS_Cpp::Transp trans_rhs);

// /////////////////////////////
// Matrix-DVector multiplication

/// v_lhs = T_M * vs_lhs (templated matrix type T_M)
template<class T_Matrix>
void V_MtV(DVector& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2);

/// vs_lhs = T_M * vs_lhs (templated matrix type T_M)
template<class T_Matrix>
void V_MtV(DVectorSlice& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2);

//@}

}

#endif // MATRIX_VECTOR_TEMPLATE_OP_H
