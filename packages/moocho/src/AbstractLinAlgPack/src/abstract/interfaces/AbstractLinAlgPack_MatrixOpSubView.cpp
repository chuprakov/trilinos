// //////////////////////////////////////////////////////////////
// MatrixWithOpSubView.cpp
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

#include <assert.h>

#include <typeinfo>
#include <stdexcept>

#include "AbstractLinAlgPack/src/MatrixWithOpSubView.h"
#include "AbstractLinAlgPack/src/MultiVectorMutable.h"
#include "AbstractLinAlgPack/src/VectorSpace.h"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/src/SpVectorClass.h"
#include "AbstractLinAlgPack/src/SpVectorView.h"
#include "AbstractLinAlgPack/src/EtaVector.h"
#include "AbstractLinAlgPack/src/LinAlgOpPack.h"
#include "ref_count_ptr.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

MatrixWithOpSubView::MatrixWithOpSubView(
	const mat_ptr_t& M_full
	,const Range1D& rng_rows
	,const Range1D& rng_cols
	,BLAS_Cpp::Transp M_trans
	)
{
	this->initialize(M_full,rng_rows,rng_cols,M_trans);
}
		
void MatrixWithOpSubView::initialize(
	const mat_ptr_t& M_full
	,const Range1D& rng_rows_in
	,const Range1D& rng_cols_in
	,BLAS_Cpp::Transp M_trans
	)
{
	namespace rcp = MemMngPack;

	if( M_full.get() ) {
		const index_type
			M_rows = M_full->rows(),
			M_cols = M_full->cols();
		const Range1D
			rng_rows = RangePack::full_range(rng_rows_in,1,M_rows),
			rng_cols = RangePack::full_range(rng_cols_in,1,M_cols);
		THROW_EXCEPTION(
			rng_rows.ubound() > M_rows, std::invalid_argument
			,"MatrixWithOpSubView::initialize(...): Error, "
			"rng_rows = ["<<rng_rows.lbound()<<","<<rng_rows.ubound()<<"] is of range of "
			"[1,M_full->rows()] = [1,"<<M_rows<<"]" );
		THROW_EXCEPTION(
			rng_cols.ubound() > M_cols, std::invalid_argument
			,"MatrixWithOpSubView::initialize(...): Error, "
			"rng_cols = ["<<rng_cols.lbound()<<","<<rng_cols.ubound()<<"] is of range of "
			"[1,M_full->cols()] = [1,"<<M_cols<<"]" );
		M_full_     = M_full;
		rng_rows_   = rng_rows;
		rng_cols_   = rng_cols;
		M_trans_    = M_trans;
		space_cols_ = ( M_trans == BLAS_Cpp::no_trans
						? M_full->space_cols().sub_space(rng_rows)->clone()
						: M_full->space_rows().sub_space(rng_cols)->clone() );
		space_rows_ = ( M_trans == BLAS_Cpp::no_trans
						? M_full->space_rows().sub_space(rng_cols)->clone()
						: M_full->space_cols().sub_space(rng_rows)->clone() );
	}
	else {
		M_full_     = rcp::null;
		rng_rows_   = Range1D::Invalid;
		rng_cols_   = Range1D::Invalid;
		M_trans_    = BLAS_Cpp::no_trans;
		space_cols_ = rcp::null;
		space_rows_ = rcp::null;
	}
}

// overridden from MatrixBase

size_type MatrixWithOpSubView::rows() const
{
	return ( M_full_.get() 
			 ? BLAS_Cpp::rows( rng_rows_.size(), rng_cols_.size(), M_trans_ )
			 : 0 );
}

size_type MatrixWithOpSubView::cols() const
{
	return ( M_full_.get() 
			 ? BLAS_Cpp::cols( rng_rows_.size(), rng_cols_.size(), M_trans_ )
			 : 0 );
}

size_type MatrixWithOpSubView::nz() const
{
	return ( M_full_.get()
			 ? ( rng_rows_.full_range() && rng_cols_.full_range()
			   ? M_full_->nz()
			   : MatrixBase::nz() )
			 : 0 );
}

// Overridden form MatrixWithOp

const VectorSpace& MatrixWithOpSubView::space_cols() const
{
	assert_initialized();
	return *space_cols_;
}

const VectorSpace& MatrixWithOpSubView::space_rows() const
{
	assert_initialized();
	return *space_rows_;
}

MatrixWithOp::mat_ptr_t
MatrixWithOpSubView::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	assert_initialized();
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

void MatrixWithOpSubView::zero_out()
{
	assert_initialized();
	if( rng_rows_.full_range() && rng_cols_.full_range() ) {
		M_full_->zero_out();
		return;
	}
	THROW_EXCEPTION(
		true, std::logic_error, "MatrixWithOpSubView::zero_out(): "
		"Error, this method can not be implemented with a sub-view" );
}

void MatrixWithOpSubView::Mt_S( value_type alpha )
{
	assert_initialized();
	if( rng_rows_.full_range() && rng_cols_.full_range() ) {
		M_full_->Mt_S(alpha);
		return;
	}
	THROW_EXCEPTION(
		true, std::logic_error, "MatrixWithOpSubView::Mt_S(alpha): "
		"Error, this method can not be implemented with a sub-view" );
}

MatrixWithOp& MatrixWithOpSubView::operator=(const MatrixWithOp& M)
{
	assert_initialized();
	assert(0); // ToDo: Implement!
	return *this;
}

std::ostream& MatrixWithOpSubView::output(std::ostream& out) const
{
	assert_initialized();
	return MatrixWithOp::output(out); // ToDo: Specialize if needed?
}

// Level-1 BLAS

// rhs matrix argument

bool MatrixWithOpSubView::Mp_StM(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	assert_initialized();
	return MatrixWithOp::Mp_StM(m_lhs,alpha,trans_rhs); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	assert_initialized();
	return MatrixWithOp::Mp_StMtP(m_lhs,alpha,M_trans,P_rhs,P_rhs_trans); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StPtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	assert_initialized();
	return MatrixWithOp::Mp_StPtM(m_lhs,alpha,P_rhs,P_rhs_trans,M_trans); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StPtMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	assert_initialized();
	return MatrixWithOp::Mp_StPtMtP(
		m_lhs,alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans); // ToDo: Specialize?
}

// lhs matrix argument

bool MatrixWithOpSubView::Mp_StM(
	value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
	assert_initialized();
	return MatrixWithOp::Mp_StM(alpha,M_rhs,trans_rhs); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StMtP(
	value_type alpha
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	assert_initialized();
	return MatrixWithOp::Mp_StMtP(alpha,M_rhs,M_trans,P_rhs,P_rhs_trans); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StPtM(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	assert_initialized();
	return MatrixWithOp::Mp_StPtM(
		alpha,P_rhs,P_rhs_trans,M_rhs,M_trans); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StPtMtP(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	assert_initialized();
	return MatrixWithOp::Mp_StPtMtP(
		alpha,P_rhs1,P_rhs1_trans,M_rhs,M_trans,P_rhs2,P_rhs2_trans); // ToDo: Specialize?
}

// Level-2 BLAS

void MatrixWithOpSubView::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans_in
	, const VectorWithOp& x, value_type b
	) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;

	assert_initialized();

	BLAS_Cpp::Transp
		M_trans_trans = ( M_trans_==no_trans ? M_trans_in : BLAS_Cpp::trans_not(M_trans_in) );

	// ToDo: Assert compatibility!

	if( rng_rows_.full_range() && rng_cols_.full_range() ) {
		AbstractLinAlgPack::Vp_StMtV(  // The matrix is just transposed
			y, a
			,*M_full_, M_trans_trans
			,x, b );
		return;
	}
	// y = b*y
	Vt_S( y, b );
	//
	// xt1                      = 0.0
	// xt3 = xt(op_op_rng_cols) = x
	// xt3                      = 0.0
	//
	// [ yt1 ]                        [ xt1 ]
	// [ yt2 ] = a * op(op(M_full)) * [ xt3 ]
	// [ yt3 ]                        [ xt3 ]
	//
	// =>
	//
	// y += yt2 = yt(op_op_rng_rows)
	//
	const Range1D
		op_op_rng_rows = ( M_trans_trans == no_trans ? rng_rows_ : rng_cols_ ),
		op_op_rng_cols = ( M_trans_trans == no_trans ? rng_cols_ : rng_rows_ );
	VectorSpace::vec_mut_ptr_t
		xt = ( M_trans_trans == no_trans ? M_full_->space_rows() : M_full_->space_cols() ).create_member(),
		yt = ( M_trans_trans == no_trans ? M_full_->space_cols() : M_full_->space_rows() ).create_member();
	*xt = 0.0;
	*xt->sub_view(op_op_rng_cols) = x;
    LinAlgOpPack::V_StMtV( yt.get(), a, *M_full_, M_trans_trans, *xt );
	LinAlgOpPack::Vp_V( y, *yt->sub_view(op_op_rng_rows) );
}

void MatrixWithOpSubView::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	assert_initialized();
	MatrixWithOp::Vp_StMtV(v_lhs,alpha,trans_rhs1,sv_rhs2,beta); // ToDo: Specialize?
}

void MatrixWithOpSubView::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& v_rhs3, value_type beta) const
{
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(
		v_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,v_rhs3,beta); // ToDo: Specialize?
}

void MatrixWithOpSubView::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(
		v_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta); // ToDo: Specialize?
}

value_type MatrixWithOpSubView::transVtMtV(
	const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorWithOp& v_rhs3) const
{
	assert_initialized();
	return MatrixWithOp::transVtMtV(v_rhs1,trans_rhs2,v_rhs3); // ToDo: Specialize?
}

value_type MatrixWithOpSubView::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	assert_initialized();
	return MatrixWithOp::transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3); // ToDo: Specialize?
}

void MatrixWithOpSubView::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert_initialized();
	MatrixWithOp::syr2k(
		M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,sym_lhs); // ToDo: Specialize?
}

// Level-3 BLAS

bool MatrixWithOpSubView::Mp_StMtM(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert_initialized();
	return MatrixWithOp::Mp_StMtM(
		m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StMtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return MatrixWithOp::Mp_StMtM(
		m_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixWithOpSubView::Mp_StMtM(
	value_type alpha
	,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
	,value_type beta )
{
	assert_initialized();
	return MatrixWithOp::Mp_StMtM(
		alpha,mvw_rhs1,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixWithOpSubView::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert_initialized();
	return MatrixWithOp::syrk(M_trans,alpha,beta,sym_lhs); // ToDo: Specialize?
}

// private

void MatrixWithOpSubView::assert_initialized() const {
	THROW_EXCEPTION(
		M_full_.get() == NULL, std::logic_error
		,"Error, the MatrixWithOpSubView object has not been initialize!" );
}

} // end namespace AbstractLinAlgPack
