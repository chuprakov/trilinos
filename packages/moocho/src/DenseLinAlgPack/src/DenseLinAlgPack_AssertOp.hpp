// //////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_AssertOp.hpp
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

#ifndef LIN_ALG_PACK_ASSERT_OP_H
#define LIN_ALG_PACK_ASSERT_OP_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

#ifdef LINALGPACK_CHECK_RHS_SIZES


/* * @name Assertion functions for linear algebra operations.
  *
  * These functions check the sizes of the linear algebra
  * expressions and throw a std::length_error if
  * the sizes do not match.  These functions
  * only perform there operations if #LINALGPACK_CHECK_RHS_SIZES#
  * is defined.
  */
// @{

/* * @name Level 1 BLAS
  */
// @{

/// v_lhs += op v_rhs
void Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size);

/// v_rhs1 op v_rhs2
void VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size);

/// op(m_lhs) += op op(m_rhs)
void Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs);

/// v_rhs1 op v_rhs2
void MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

//		end Level 1 BLAS
// @}

/* *  @name Level 2 BLAS
  */
// @{

/// op(m_rhs1) * v_rhs2
void MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
	, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size);

/// v_lhs += op(m_rhs1) * v_rhs2
void Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
	, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size);

//		end Level 2 BLAS
// @}

/* *  @name Level 3 BLAS
  */
// @{

/// op(m_lhs) += op(m_rhs1)
void MtM_assert_sizes(
	  size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

/// op(m_lhs) += op(m_rhs1) * op(m_rhs2)
void Mp_MtM_assert_sizes(
	  size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

//		end Level 3 BLAS
// @}


// @}

#else

// inline definitions that do nothing

inline
void Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size)
{}

inline
void VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size)
{}

inline
void Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs)
{}

inline
void MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

inline
void MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
	, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{}

inline
void Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
	, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{}

inline
void MtM_assert_sizes(
	  size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

inline
void Mp_MtM_assert_sizes(
	  size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

#endif

} // end namespace DenseLinAlgPack

#endif	// LIN_ALG_PACK_ASSERT_OP_H
