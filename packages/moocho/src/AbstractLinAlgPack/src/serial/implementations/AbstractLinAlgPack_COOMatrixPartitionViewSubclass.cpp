// //////////////////////////////////////////////////////////////////////////////////
// COOMatrixPartitionViewSubclass.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include "LinAlgPack/include/LinAlgOpPack.h"
#include "../include/COOMatrixPartitionViewSubclass.h"
#include "../include/SparseVectorSliceOp.h"
#include "../include/SparseElement.h"
#include "../include/COOMPartitionOp.h"
#include "LinAlgPack/include/GenMatrixOp.h"

namespace LinAlgOpPack {

using SparseLinAlgPack::Vp_StV;
using SparseLinAlgPack::Vp_StMtV;
using SparseLinAlgPack::Mp_StM;
using SparseLinAlgPack::Mp_StMtM;

}	// end namespace LinAlgOpPack

namespace SparseLinAlgPack {

size_type COOMatrixPartitionViewSubclass::rows() const {
	return trans_ == BLAS_Cpp::no_trans ? m().rows() : m().cols();
}

size_type COOMatrixPartitionViewSubclass::cols() const {
	return trans_ == BLAS_Cpp::no_trans ? m().cols() : m().rows();
}

MatrixWithOp& COOMatrixPartitionViewSubclass::operator=(const MatrixWithOp& m) {
	if(&m == this) return *this;	// assignment to self
	const COOMatrixPartitionViewSubclass *p_m = dynamic_cast<const COOMatrixPartitionViewSubclass*>(&m);
	if(p_m) {
		throw std::invalid_argument("COOMatrixPartitionViewSubclass::operator=(const MatrixWithOp& m)"
			" :  There is not an assignment operator defined for COOMatrixWithPartitionedView::partition_type"
			".   Only assignment to self can be handeled" );
	}
	else {
		throw std::invalid_argument("COOMatrixPartitionViewSubclass::operator=(const MatrixWithOp& m)"
			" : The concrete type of m is not a subclass of COOMatrixPartitionViewSubclass as expected" );
	}
	return *this;
}

// Level-1 BLAS

void COOMatrixPartitionViewSubclass::Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	SparseLinAlgPack::Mp_StM(gms_lhs,alpha,m(),op(trans_rhs));
}

// Level-2 BLAS

void COOMatrixPartitionViewSubclass::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta) const
{
	SparseLinAlgPack::Vp_StMtV(vs_lhs, alpha, m(), op(trans_rhs1), vs_rhs2, beta);
}

void COOMatrixPartitionViewSubclass::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	Vector v_rhs2;
	LinAlgOpPack::assign(&v_rhs2,sv_rhs2);
	SparseLinAlgPack::Vp_StMtV(vs_lhs, alpha, m(), op(trans_rhs1), v_rhs2(), beta);
}

value_type COOMatrixPartitionViewSubclass::transVtMtV(const VectorSlice& vs_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const VectorSlice& vs_rhs3) const
{
	Vector tmp;
	LinAlgOpPack::V_MtV(&tmp,m(),op(trans_rhs2),vs_rhs3);
	return LinAlgPack::dot(vs_rhs1,tmp());
}

value_type COOMatrixPartitionViewSubclass::transVtMtV(const SpVectorSlice& sv_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
	Vector v_rhs3;
	LinAlgOpPack::assign(&v_rhs3,sv_rhs3);
	Vector tmp;
	LinAlgOpPack::V_MtV(&tmp,m(),op(trans_rhs2),v_rhs3());
	return dot(sv_rhs1,tmp());
}

// Level-3 BLAS

void COOMatrixPartitionViewSubclass::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	SparseLinAlgPack::Mp_StMtM(gms_lhs, alpha, m(), op(trans_rhs1), gms_rhs2, trans_rhs2, beta);
}

void COOMatrixPartitionViewSubclass::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	SparseLinAlgPack::Mp_StMtM(gms_lhs, alpha, gms_rhs1, trans_rhs1, m(), op(trans_rhs2), beta);
}

}	// end namespace SparseLinAlgPack