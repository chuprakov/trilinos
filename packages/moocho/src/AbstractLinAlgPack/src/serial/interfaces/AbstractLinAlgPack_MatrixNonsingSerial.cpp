// //////////////////////////////////////////////////////////////////////////////////
//  MatrixNonsingularSerial.cpp
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

// ToDo: 3/6/00: Provide default implementations for these
// operations.

#include <assert.h>

#include "SparseLinAlgPack/src/MatrixNonsingularSerial.hpp"
#include "SparseLinAlgPack/src/MatrixWithOpSerial.hpp"
#include "SparseLinAlgPack/src/VectorDenseEncap.hpp"
#include "SparseLinAlgPack/src/VectorWithOpGetSparse.hpp"
#include "SparseLinAlgPack/src/MatrixWithOpGetGMSMutable.hpp"
#include "SparseLinAlgPack/src/MatrixWithOpGetGMSTri.hpp"
#include "SparseLinAlgPack/src/MatrixSymWithOpGetGMSSymMutable.hpp"
#include "SparseLinAlgPack/src/SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/SpVectorClass.hpp"
#include "LinAlgPack/src/GenMatrixClass.hpp"
#include "LinAlgPack/src/VectorClass.hpp"
#include "LinAlgPack/src/LinAlgOpPack.hpp"
#include "LinAlgPack/src/LinAlgPackAssertOp.hpp"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Mp_StM;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace SparseLinAlgPack {

//  Level-2 BLAS

void MatrixNonsingularSerial::V_InvMtV(
	Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1,const VectorSlice& vs_rhs2
	) const
{
	const size_type n = rows();
	LinAlgPack::MtV_assert_sizes( n, n, trans_rhs1, vs_rhs2.dim() );
	v_lhs->resize(n);
	this->V_InvMtV( &(*v_lhs)(), trans_rhs1, vs_rhs2 );
}

void MatrixNonsingularSerial::V_InvMtV(
	Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
	) const
{
	const size_type n = rows();
	LinAlgPack::MtV_assert_sizes( n, n, trans_rhs1, sv_rhs2.dim() );
	v_lhs->resize(n);
	Vector v_rhs2;
	LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
	this->V_InvMtV( &(*v_lhs)(), trans_rhs1, v_rhs2() );
}

void MatrixNonsingularSerial::V_InvMtV(
	VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
	) const
{
	const size_type n = rows();
	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->dim(), n, n, trans_rhs1, sv_rhs2.dim() );
	Vector v_rhs2;
	LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
	this->V_InvMtV( vs_lhs, trans_rhs1, v_rhs2() );
}

value_type MatrixNonsingularSerial::transVtInvMtV(
	const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2, const VectorSlice& vs_rhs3
	) const
{
	const size_type n = rows();
	LinAlgPack::Vp_MtV_assert_sizes( vs_rhs1.dim(), n, n, trans_rhs2, vs_rhs3.dim() );
	Vector tmp;
	this->V_InvMtV( &tmp, trans_rhs2, vs_rhs3 );
	return LinAlgPack::dot( vs_rhs1, tmp() );
}

value_type MatrixNonsingularSerial::transVtInvMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3
	) const
{
	const size_type n = rows();
	LinAlgPack::Vp_MtV_assert_sizes( sv_rhs1.dim(), n, n, trans_rhs2, sv_rhs3.dim() );
	Vector tmp;
	this->V_InvMtV( &tmp, trans_rhs2, sv_rhs3 );
	return SparseLinAlgPack::dot( sv_rhs1, tmp() );
}

// Level-3 BLAS

void MatrixNonsingularSerial::M_StInvMtM(
	GenMatrix* C, value_type a
	,BLAS_Cpp::Transp A_trans
	,const GenMatrixSlice& B, BLAS_Cpp::Transp B_trans
	) const
{
	LinAlgPack::MtM_assert_sizes( rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	C->resize(
		  BLAS_Cpp::rows( rows(), cols(), A_trans )
		, BLAS_Cpp::cols( B.rows(), B.cols(), B_trans )
		);
	this->M_StInvMtM( &(*C)(), a, A_trans, B, B_trans );
}

void MatrixNonsingularSerial::M_StInvMtM(
	GenMatrixSlice* C, value_type a
	,BLAS_Cpp::Transp A_trans
	,const GenMatrixSlice& B, BLAS_Cpp::Transp B_trans
	) const
{
	LinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	//
	// C = a * inv(op(A)) * op(B)
	//
	// C.col(j) = a * inv(op(A)) * op(B).col(j)
	//

	for( size_type j = 1; j <= C->cols(); ++j )
		SparseLinAlgPack::V_InvMtV( &C->col(j), *this, A_trans
			, LinAlgPack::col( B, B_trans, j ) );
	if( a != 1.0 )
		LinAlgOpPack::Mt_S( C, a );
}

void MatrixNonsingularSerial::M_StMtInvM(
	GenMatrix* gm_lhs, value_type alpha
	,const GenMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	assert(0);	// ToDo: Implement this!
}

void MatrixNonsingularSerial::M_StMtInvM(
	GenMatrixSlice* gms_lhs, value_type alpha
	,const GenMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	assert(0);	// ToDo: Implement this!
}

void MatrixNonsingularSerial::M_StInvMtM(
	GenMatrix* C, value_type a
	,BLAS_Cpp::Transp A_trans
	,const MatrixWithOpSerial& B, BLAS_Cpp::Transp B_trans
	) const
{
	LinAlgPack::MtM_assert_sizes( rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	C->resize(
		  BLAS_Cpp::rows( rows(), cols(), A_trans )
		, BLAS_Cpp::cols( B.rows(), B.cols(), B_trans )
		);
	SparseLinAlgPack::M_StInvMtM( &(*C)(), a, *this, A_trans, B, B_trans );
}

void MatrixNonsingularSerial::M_StInvMtM(
	GenMatrixSlice* C, value_type a
	,BLAS_Cpp::Transp A_trans
	,const MatrixWithOpSerial& B, BLAS_Cpp::Transp B_trans
	) const
{
	using LinAlgOpPack::assign;
	LinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	GenMatrix B_dense;
	assign( &B_dense, B, BLAS_Cpp::no_trans );
	SparseLinAlgPack::M_StInvMtM( C, a, *this, A_trans, B_dense(), B_trans );
}

void MatrixNonsingularSerial::M_StMtInvM(
	GenMatrix* gm_lhs, value_type alpha
	,const MatrixWithOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	assert(0);	// ToDo: Implement this!
}

void MatrixNonsingularSerial::M_StMtInvM(
	GenMatrixSlice* gms_lhs, value_type alpha
	,const MatrixWithOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	assert(0);	// ToDo: Implement this!
}

// Overridden from MatrixNonsingular

void MatrixNonsingularSerial::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const VectorWithOp& v_rhs2) const
{
	VectorDenseMutableEncap       vs_lhs(*v_lhs);
	const VectorWithOpGetSparse   *sv_rhs2 = dynamic_cast<const VectorWithOpGetSparse*>(&v_rhs2);
	if(sv_rhs2)
		this->V_InvMtV( &vs_lhs(), trans_rhs1, VectorSparseEncap(*sv_rhs2)() );
	VectorDenseEncap              vs_rhs2(v_rhs2);
	this->V_InvMtV( &vs_lhs(), trans_rhs1, vs_rhs2() );	
}

void MatrixNonsingularSerial::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const SpVectorSlice& sv_rhs2) const
{
	this->V_InvMtV( &VectorDenseMutableEncap(*v_lhs)(), trans_rhs1, sv_rhs2 );
}

value_type MatrixNonsingularSerial::transVtInvMtV(
	const VectorWithOp& v_rhs1
	,BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3) const
{
	VectorDenseEncap              vs_rhs1(v_rhs1);
	VectorDenseEncap              vs_rhs3(v_rhs3);
	return this->transVtInvMtV(vs_rhs1(),trans_rhs2,vs_rhs3());
}

void MatrixNonsingularSerial::M_StInvMtM(
	MatrixWithOp* m_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	MatrixDenseMutableEncap
		gms_lhs(m_lhs);      // Warning!  This may throw an exception!
	if(const MatrixWithOpGetGMS* mwo_gms_rhs2 = dynamic_cast<const MatrixWithOpGetGMS*>(&mwo_rhs2)) {
		this->M_StInvMtM(&gms_lhs(),alpha,trans_rhs1,MatrixDenseEncap(*mwo_gms_rhs2)(),trans_rhs2);
		return;
	}
	this->M_StInvMtM(&gms_lhs(),alpha,trans_rhs1,dyn_cast<const MatrixWithOpSerial>(mwo_rhs2),trans_rhs2);
}

void MatrixNonsingularSerial::M_StMtInvM(
	MatrixWithOp* m_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	MatrixDenseMutableEncap
		gms_lhs(m_lhs);      // Warning!  This may throw an exception!
	if(const MatrixWithOpGetGMS* mwo_gms_rhs1 = dynamic_cast<const MatrixWithOpGetGMS*>(&mwo_rhs1)) {
		this->M_StMtInvM(&gms_lhs(),alpha,MatrixDenseEncap(*mwo_gms_rhs1)(),trans_rhs1,trans_rhs2);
		return;
	}
	this->M_StMtInvM(&gms_lhs(),alpha,dyn_cast<const MatrixWithOpSerial>(mwo_rhs1),trans_rhs1,trans_rhs2);
}

}	// end namespace SparseLinAlgPack
