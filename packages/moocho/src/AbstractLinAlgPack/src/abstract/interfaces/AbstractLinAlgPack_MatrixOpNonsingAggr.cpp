// /////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingularAggr.cpp
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

#include "AbstractLinAlgPack/include/MatrixWithOpNonsingularAggr.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace AbstractLinAlgPack {

// Constructors / initializers

MatrixWithOpNonsingularAggr::MatrixWithOpNonsingularAggr()
{} // Nothing to explicitly initialize

MatrixWithOpNonsingularAggr::MatrixWithOpNonsingularAggr(
	const mwo_ptr_t       &mwo
	,BLAS_Cpp::Transp     mwo_trans
	,const mns_ptr_t      &mns
	,BLAS_Cpp::Transp     mns_trans
	)
{
	this->initialize(mwo,mwo_trans,mns,mns_trans);
}

void MatrixWithOpNonsingularAggr::initialize(
	const mwo_ptr_t       &mwo
	,BLAS_Cpp::Transp     mwo_trans
	,const mns_ptr_t      &mns
	,BLAS_Cpp::Transp     mns_trans
	)
{
#ifdef _DEBUG
	const size_type
		mwo_rows = mwo_->rows(),
		mwo_cols = mwo_->cols(),
		mns_rows = mns_->rows(),
		mns_cols = mns_->cols();
	THROW_EXCEPTION(
		mwo.get() == NULL, std::invalid_argument
		,"MatrixWithOpNonsingularAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mns.get() == NULL, std::invalid_argument
		,"MatrixWithOpNonsingularAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mwo_rows == (mwo_trans == mns_trans ? mns_rows : mns_cols), std::invalid_argument
		,"MatrixWithOpNonsingularAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mwo_cols == (mwo_trans == mns_trans ? mns_cols : mns_rows), std::invalid_argument
		,"MatrixWithOpNonsingularAggr::initialize(...): Error!" );
#endif
	mwo_       = mwo;
	mwo_trans_ = mwo_trans;
	mns_       = mns;
	mns_trans_ = mns_trans;
	THROW_EXCEPTION(
		mwo_trans_ == BLAS_Cpp::trans || mns_trans_ == BLAS_Cpp::trans, std::logic_error
		,"MatrixWithOpNonsingularAggr::initialize(...): Error, "
		"Can't handle transposed matrices yet!" );
}

void MatrixWithOpNonsingularAggr::set_uninitialized()
{
	namespace rcp = MemMngPack;
	mwo_       = rcp::null;
	mwo_trans_ = BLAS_Cpp::no_trans;
	mns_       = rcp::null;
	mns_trans_ = BLAS_Cpp::no_trans;
}

// Overridden from MatrixBase

size_type MatrixWithOpNonsingularAggr::rows() const
{
	return mwo_->rows(); // square matrix!
}

size_type MatrixWithOpNonsingularAggr::cols() const
{
	return mwo_->cols(); // square matrix!
}

size_type MatrixWithOpNonsingularAggr::nz() const
{
	return mwo_->nz();
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixWithOpNonsingularAggr::space_cols() const
{
	return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_cols() : mwo_->space_rows();
}

const VectorSpace& MatrixWithOpNonsingularAggr::space_rows() const
{
	return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_rows() : mwo_->space_cols();
}

MatrixWithOp::mat_ptr_t
MatrixWithOpNonsingularAggr::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	return MatrixWithOp::sub_view(row_rng,col_rng); // ToDo: Speicalize!
}

MatrixWithOp& MatrixWithOpNonsingularAggr::operator=(const MatrixWithOp& M)
{
	using DynamicCastHelperPack::dyn_cast;
	const MatrixWithOpNonsingularAggr
		Mp = dyn_cast<const MatrixWithOpNonsingularAggr>(M);
	if( this == &Mp )
		return *this; // Assignment to self
	// Shallow copy is okay as long as client is careful!
	mwo_       = Mp.mwo_;
	mwo_trans_ = Mp.mwo_trans_;
	mns_       = Mp.mns_;
	mns_trans_ = Mp.mns_trans_;
}

std::ostream& MatrixWithOpNonsingularAggr::output(std::ostream& out) const
{
	out << "Aggregate nonsingular matrix:\n";
	out << "mwo =\n" << *mwo_;
	out << "mns = ???\n";
	return out;
}

// ToDo: Update below code for transposed matrices!

bool MatrixWithOpNonsingularAggr::Mp_StM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	return mwo_->Mp_StM(mwo_lhs,alpha,trans_rhs);
}

bool MatrixWithOpNonsingularAggr::Mp_StMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	return mwo_->Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans);
}

bool MatrixWithOpNonsingularAggr::Mp_StPtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	return mwo_->Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans);
}

bool MatrixWithOpNonsingularAggr::Mp_StPtMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	return mwo_->Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans);
}

void MatrixWithOpNonsingularAggr::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorWithOp& x, value_type b) const
{
	mwo_->Vp_StMtV(y,a,M_trans,x,b);
}

void MatrixWithOpNonsingularAggr::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	mwo_->Vp_StMtV(y,a,M_trans,x,b);
}

void MatrixWithOpNonsingularAggr::Vp_StPtMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& v_rhs3, value_type beta) const
{
	mwo_->Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,v_rhs3,beta);
}

void MatrixWithOpNonsingularAggr::Vp_StPtMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	mwo_->Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta);
}

value_type MatrixWithOpNonsingularAggr::transVtMtV(
	const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorWithOp& v_rhs3) const
{
	return mwo_->transVtMtV(v_rhs1,trans_rhs2,v_rhs3);
}

value_type MatrixWithOpNonsingularAggr::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	return mwo_->transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

void MatrixWithOpNonsingularAggr::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymWithOp* symwo_lhs ) const
{
	mwo_->syr2k(M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs);
}

bool MatrixWithOpNonsingularAggr::Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return mwo_->Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

bool MatrixWithOpNonsingularAggr::Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return mwo_->Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta);
}

void MatrixWithOpNonsingularAggr::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
		mwo_->syrk(M_trans,alpha,beta,sym_lhs);
}

// Overridden from MatrixNonsingular */

// ToDo: Update below code for transposed matrices!

void MatrixWithOpNonsingularAggr::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const VectorWithOp& v_rhs2) const
{
	mns_->V_InvMtV(v_lhs,trans_rhs1,v_rhs2);
}

void MatrixWithOpNonsingularAggr::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	mns_->V_InvMtV(v_lhs,trans_rhs1,sv_rhs2);
}

value_type MatrixWithOpNonsingularAggr::transVtInvMtV(
	const VectorWithOp& v_rhs1
	,BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3) const
{
	return mns_->transVtInvMtV(v_rhs1,trans_rhs2,v_rhs3);
}

value_type MatrixWithOpNonsingularAggr::transVtInvMtV(
	const SpVectorSlice& sv_rhs1
	,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
	return mns_->transVtInvMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

void MatrixWithOpNonsingularAggr::M_StInvMtM(
	MatrixWithOp* m_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	) const
{
	mns_->M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

void MatrixWithOpNonsingularAggr::M_StMtInvM(
	MatrixWithOp* m_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	mns_->M_StMtInvM(m_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

} // end namespace AbstractLinAlgPack
