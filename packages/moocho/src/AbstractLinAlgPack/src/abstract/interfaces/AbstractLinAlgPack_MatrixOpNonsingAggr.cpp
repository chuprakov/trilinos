// /////////////////////////////////////////////////////////////////////
// MatrixOpNonsingAggr.cpp
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

#include "AbstractLinAlgPack/src/MatrixOpNonsingAggr.hpp"
#include "AbstractLinAlgPack/src/MatrixOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "ThrowException.hpp"
#include "dynamic_cast_verbose.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

MatrixOpNonsingAggr::MatrixOpNonsingAggr()
{} // Nothing to explicitly initialize

MatrixOpNonsingAggr::MatrixOpNonsingAggr(
	const mwo_ptr_t       &mwo
	,BLAS_Cpp::Transp     mwo_trans
	,const mns_ptr_t      &mns
	,BLAS_Cpp::Transp     mns_trans
	)
{
	this->initialize(mwo,mwo_trans,mns,mns_trans);
}

void MatrixOpNonsingAggr::initialize(
	const mwo_ptr_t       &mwo
	,BLAS_Cpp::Transp     mwo_trans
	,const mns_ptr_t      &mns
	,BLAS_Cpp::Transp     mns_trans
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		mwo.get() == NULL, std::invalid_argument
		,"MatrixOpNonsingAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mns.get() == NULL, std::invalid_argument
		,"MatrixOpNonsingAggr::initialize(...): Error!" );
	const size_type
		mwo_rows = mwo->rows(),
		mwo_cols = mwo->cols(),
		mns_rows = mns->rows(),
		mns_cols = mns->cols();
	THROW_EXCEPTION(
		mwo_rows != mwo_cols, std::invalid_argument
		,"MatrixOpNonsingAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mns_rows != mns_cols, std::invalid_argument
		,"MatrixOpNonsingAggr::initialize(...): Error!" );
	THROW_EXCEPTION(
		mwo_rows != mns_rows, std::invalid_argument
		,"MatrixOpNonsingAggr::initialize(...): Error!" );
#endif
	mwo_       = mwo;
	mwo_trans_ = mwo_trans;
	mns_       = mns;
	mns_trans_ = mns_trans;
}

void MatrixOpNonsingAggr::set_uninitialized()
{
	namespace rcp = MemMngPack;
	mwo_       = rcp::null;
	mwo_trans_ = BLAS_Cpp::no_trans;
	mns_       = rcp::null;
	mns_trans_ = BLAS_Cpp::no_trans;
}

// Overridden from MatrixBase

size_type MatrixOpNonsingAggr::rows() const
{
	return mwo_.get() ? mwo_->rows() : 0; // square matrix!
}

size_type MatrixOpNonsingAggr::cols() const
{
	return mwo_.get() ? mwo_->rows() : 0; // square matrix!
}

size_type MatrixOpNonsingAggr::nz() const
{
	return mwo_.get() ? mwo_->nz() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixOpNonsingAggr::space_cols() const
{
	return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_cols() : mwo_->space_rows();
}

const VectorSpace& MatrixOpNonsingAggr::space_rows() const
{
	return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_rows() : mwo_->space_cols();
}

MatrixOp::mat_ptr_t
MatrixOpNonsingAggr::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	return MatrixOp::sub_view(row_rng,col_rng); // ToDo: Speicalize!
}

MatrixOp& MatrixOpNonsingAggr::operator=(const MatrixOp& M)
{
	using DynamicCastHelperPack::dyn_cast;
	const MatrixOpNonsingAggr
		Mp = dyn_cast<const MatrixOpNonsingAggr>(M);
	if( this == &Mp )
		return *this; // Assignment to self
	// Shallow copy is okay as long as client is careful!
	mwo_       = Mp.mwo_;
	mwo_trans_ = Mp.mwo_trans_;
	mns_       = Mp.mns_;
	mns_trans_ = Mp.mns_trans_;
	return *this;
}

std::ostream& MatrixOpNonsingAggr::output(std::ostream& out) const
{
	out << "Aggregate nonsingular matrix:\n";
	out << (mwo_trans_ == BLAS_Cpp::no_trans ? "mwo" : "mwo\'") << " =\n" << *mwo_;
	out << (mns_trans_ == BLAS_Cpp::no_trans ? "mns" : "mns\'") << " = ???\n";
	return out;
}

bool MatrixOpNonsingAggr::Mp_StM(
	MatrixOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	return mwo_->Mp_StM(mwo_lhs,alpha,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs));
}

bool MatrixOpNonsingAggr::Mp_StMtP(
	MatrixOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	return mwo_->Mp_StMtP(
		mwo_lhs,alpha,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),P_rhs,P_rhs_trans);
}

bool MatrixOpNonsingAggr::Mp_StPtM(
	MatrixOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	return mwo_->Mp_StPtM(
		mwo_lhs,alpha,P_rhs,P_rhs_trans,BLAS_Cpp::trans_trans(mwo_trans_,M_trans));
}

bool MatrixOpNonsingAggr::Mp_StPtMtP(
	MatrixOp* mwo_lhs, value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	return mwo_->Mp_StPtMtP(
		mwo_lhs,alpha,P_rhs1,P_rhs1_trans,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),P_rhs2,P_rhs2_trans);
}

void MatrixOpNonsingAggr::Vp_StMtV(
	VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const Vector& x, value_type b) const
{
	mwo_->Vp_StMtV(y,a,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),x,b);
}

void MatrixOpNonsingAggr::Vp_StMtV(
	VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	mwo_->Vp_StMtV(y,a,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),x,b);
}

void MatrixOpNonsingAggr::Vp_StPtMtV(
	VectorMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const Vector& v_rhs3, value_type beta) const
{
	mwo_->Vp_StPtMtV(
		vs_lhs,alpha,P_rhs1,P_rhs1_trans,BLAS_Cpp::trans_trans(mwo_trans_,M_rhs2_trans),v_rhs3,beta);
}

void MatrixOpNonsingAggr::Vp_StPtMtV(
	VectorMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	mwo_->Vp_StPtMtV(
		vs_lhs,alpha,P_rhs1,P_rhs1_trans,BLAS_Cpp::trans_trans(mwo_trans_,M_rhs2_trans),sv_rhs3,beta);
}

value_type MatrixOpNonsingAggr::transVtMtV(
	const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const Vector& v_rhs3) const
{
	return mwo_->transVtMtV(v_rhs1,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),v_rhs3);
}

value_type MatrixOpNonsingAggr::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	return mwo_->transVtMtV(sv_rhs1,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),sv_rhs3);
}

void MatrixOpNonsingAggr::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymOp* symwo_lhs ) const
{
	mwo_->syr2k(BLAS_Cpp::trans_trans(mwo_trans_,M_trans),alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs);
}

bool MatrixOpNonsingAggr::Mp_StMtM(
	MatrixOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return mwo_->Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),beta);
}

bool MatrixOpNonsingAggr::Mp_StMtM(
	MatrixOp* mwo_lhs, value_type alpha
	, const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return mwo_->Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),beta);
}

bool MatrixOpNonsingAggr::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymOp* sym_lhs ) const
{
		return mwo_->syrk(BLAS_Cpp::trans_trans(mwo_trans_,M_trans),alpha,beta,sym_lhs);
}

// Overridden from MatrixNonsing */

void MatrixOpNonsingAggr::V_InvMtV(
	VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const Vector& v_rhs2) const
{
	mns_->V_InvMtV(v_lhs,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),v_rhs2);
}

void MatrixOpNonsingAggr::V_InvMtV(
	VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	mns_->V_InvMtV(v_lhs,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),sv_rhs2);
}

value_type MatrixOpNonsingAggr::transVtInvMtV(
	const Vector& v_rhs1
	,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3) const
{
	return mns_->transVtInvMtV(v_rhs1,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs2),v_rhs3);
}

value_type MatrixOpNonsingAggr::transVtInvMtV(
	const SpVectorSlice& sv_rhs1
	,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
	return mns_->transVtInvMtV(sv_rhs1,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs2),sv_rhs3);
}

void MatrixOpNonsingAggr::M_StInvMtM(
	MatrixOp* m_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	) const
{
	mns_->M_StInvMtM(m_lhs,alpha,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),mwo_rhs2,trans_rhs2);
}

void MatrixOpNonsingAggr::M_StMtInvM(
	MatrixOp* m_lhs, value_type alpha
	,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	mns_->M_StMtInvM(m_lhs,alpha,mwo_rhs1,trans_rhs1,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1));
}

} // end namespace AbstractLinAlgPack
