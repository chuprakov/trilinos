// ///////////////////////////////////////////////////////////////////////////
// MatrixVectorTemplateOpDef.h
//
// Definitions of template functions declared in MatrixVectorTemplateOp.h.

#ifndef MATRIX_VECTOR_TEMPLATE_OP_DEF_H
#define MATRIX_VECTOR_TEMPLATE_OP_DEF_H

#include "MatrixVectorTemplateOp.h"
#include "LinAlgPack/include/GenMatrixClass.h"

// ///////////////////////////////////
// Matrix assignment

namespace {
// Typedef for vector returning functions (row or col but not diagonal)
typedef SparseLinAlgPack::VectorSlice (SparseLinAlgPack::GenMatrixSlice::*Pvec_func)
	(SparseLinAlgPack::GenMatrixSlice::size_type);
// Implement sparse matrix, dense matrix assignment.  Sizes are not checked.
template<class T_Matrix>
void imp_assign(SparseLinAlgPack::GenMatrixSlice& gms_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	// If trans, then copy col into row, otherwise copy col into col.
	Pvec_func vec_func;
	if(trans_rhs == BLAS_Cpp::no_trans)		vec_func = &SparseLinAlgPack::GenMatrixSlice::col;
	else									vec_func = &SparseLinAlgPack::GenMatrixSlice::row;
	for(int k = 1; k <= gm_rhs.cols(); ++k)
		SparseLinAlgPack::assign((gms_lhs.*vec_func)(k), gm_rhs.col(k));
}
} // end namespace

// Definitions of template functions for matrix-matrix assignment

/// gm_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void SparseLinAlgPack::assign(GenMatrix& gm_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	LinAlgPack::resize_gm_lhs(gm_lhs,gm_rhs.rows(),gm_rhs.cols(),trans_rhs);
	GenMatrixSlice gms_lhs = gm_lhs;
	imp_assign(gms_lhs,gm_rhs,trans_rhs);
}

/// gms_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void SparseLinAlgPack::assign(GenMatrixSlice& gms_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	LinAlgPack::assert_gms_lhs(gms_lhs,gm_rhs.rows(),gm_rhs.cols(),trans_rhs);
	imp_assign(gms_lhs,gm_rhs,trans_rhs);
}

// //////////////////////////////////
// Matrix-Vector multiplication

namespace {
// Throw an exeption of the rhs arguments do not match
template<class T_Matrix>
void imp_assert_V_MtV_rhs_sizes(const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const SparseLinAlgPack::VectorSlice& vs_rhs2)
{
	typename T_Matrix::size_type
		cols = (trans_rhs1 == BLAS_Cpp::no_trans) ? gm_rhs1.cols() : gm_rhs1.rows();

	if(cols != vs_rhs2.size())
		throw std::length_error("V_MtV: The sizes of the rhs expression do not match");
}

// Implementation of matrix-vector multiply (no transpose).  Sizes are not checked.
template<class T_Matrix>
void imp_V_MtV_no_trans(SparseLinAlgPack::VectorSlice& vs_lhs, const T_Matrix& gm_rhs1
	, const SparseLinAlgPack::VectorSlice& vs_rhs2)
{
	typedef typename T_Matrix::size_type size_type;
	size_type rows = gm_rhs1.rows();
	SparseLinAlgPack::VectorSlice::iterator itr_v_lhs = vs_lhs.begin();
	for(size_type i = 1; i <= rows; ++i)
		*itr_v_lhs++ = SparseLinAlgPack::dot(vs_rhs2,gm_rhs1.row(i));
}
// Implementation of matrix-vector multiply (transpose).  Sizes are not checked.
template<class T_Matrix>
void imp_V_MtV_trans(SparseLinAlgPack::VectorSlice& vs_lhs, const T_Matrix& gm_rhs1
	, const SparseLinAlgPack::VectorSlice& vs_rhs2)
{
	typedef typename T_Matrix::size_type size_type;
	size_type cols = gm_rhs1.cols();
	SparseLinAlgPack::VectorSlice::iterator itr_v_lhs = vs_lhs.begin();
	for(size_type j = 1; j <= cols; ++j)
		*itr_v_lhs++ = SparseLinAlgPack::dot(vs_rhs2,gm_rhs1.col(j));
}

} // end namespace

// Definitions of template functions for matrix-vector multiplication

template<class T_Matrix>
void SparseLinAlgPack::V_MtV(Vector& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	imp_assert_V_MtV_rhs_sizes(gm_rhs1,trans_rhs1,vs_rhs2);
	v_lhs.resize( (trans_rhs1==BLAS_Cpp::no_trans) ? gm_rhs1.rows() : gm_rhs1.cols() );
	VectorSlice vs_lhs = v_lhs;
	if(trans_rhs1 == BLAS_Cpp::no_trans)
		imp_V_MtV_no_trans(vs_lhs,gm_rhs1,vs_rhs2);
	else
		imp_V_MtV_trans(vs_lhs,gm_rhs1,vs_rhs2);
}

template<class T_Matrix>
void SparseLinAlgPack::V_MtV(VectorSlice& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	imp_assert_V_MtV_rhs_sizes(gm_rhs1,trans_rhs1,vs_rhs2);
	LinAlgPack::assert_resize_vs_lhs(v_lhs, (trans_rhs1==BLAS_Cpp::no_trans) ? gm_rhs1.rows() : gm_rhs1.cols());
	if(trans_rhs1 == BLAS_Cpp::no_trans)
		imp_V_MtV_no_trans(v_lhs,gm_rhs1,vs_rhs2);
	else
		imp_V_MtV_trans(v_lhs,gm_rhs1,vs_rhs2);
}

#endif // MATRIX_VECTOR_TEMPLATE_OP_DEF_H
