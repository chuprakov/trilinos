// /////////////////////////////////////////////////////////////////////////////
// MultiVectorMutableCols.cpp
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

#include "SparseLinAlgPack/include/MultiVectorMutableCols.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpGetGMSSymMutable.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

// Constructors/Initializers

MultiVectorMutableCols::MultiVectorMutableCols()
{}

MultiVectorMutableCols::MultiVectorMutableCols(
	const  MemMngPack::ref_count_ptr<const VectorSpace>   &space_cols
	,const  MemMngPack::ref_count_ptr<const VectorSpace>  &space_rows
	,MemMngPack::ref_count_ptr<VectorWithOpMutable>       col_vecs[]
	)
{
	this->initialize(space_cols,space_rows,col_vecs);
}
	
void MultiVectorMutableCols::initialize(
	const  MemMngPack::ref_count_ptr<const VectorSpace>   &space_cols
	,const  MemMngPack::ref_count_ptr<const VectorSpace>  &space_rows
	,MemMngPack::ref_count_ptr<VectorWithOpMutable>       col_vecs[]
	)
{
	space_cols_ = space_cols;
	space_rows_ = space_rows;
	const size_type num_cols = space_rows->dim();
	col_vecs_.resize(num_cols);
	if(col_vecs) {
		for( size_type j = 1; j <= num_cols; ++j )
			col_vecs_[j-1] = col_vecs[j-1];
	}
	else {
		for( size_type j = 1; j <= num_cols; ++j )
			col_vecs_[j-1] = space_cols->create_member();
	}
}

void MultiVectorMutableCols::set_uninitialized()
{
	col_vecs_.resize(0);
	space_cols_ = MemMngPack::null;
	space_rows_ = MemMngPack::null;
}

// Overridden from MatrixBase

size_type MultiVectorMutableCols::rows() const
{
	return space_cols_.get() ? space_cols_->dim() : 0;
}

size_type MultiVectorMutableCols::cols() const
{
	return space_rows_.get() ? space_rows_->dim() : 0;
}

// Overridden from MatrixWithOp

const VectorSpace& MultiVectorMutableCols::space_cols() const
{
	return *space_cols_;
}

const VectorSpace& MultiVectorMutableCols::space_rows() const
{
	return *space_rows_;
}

void MultiVectorMutableCols::zero_out()
{
	for( size_type j = 1; j <= col_vecs_.size(); ++j )
		col_vecs_[j-1]->zero();
}

void MultiVectorMutableCols::Mt_S( value_type alpha )
{
	for( size_type j = 1; j <= col_vecs_.size(); ++j )
		LinAlgOpPack::Vt_S(col_vecs_[j-1].get(),alpha);
}

MatrixWithOp&
MultiVectorMutableCols::operator=(const MatrixWithOp& mwo_rhs)
{
	const MultiVector
		*mv_rhs = dynamic_cast<const MultiVector*>(&mwo_rhs);
	if(!mv_rhs)
		return MatrixWithOp::operator=(mwo_rhs);
	for( size_type j = 1; j <= col_vecs_.size(); ++j )
		*col_vecs_[j-1] = *mv_rhs->col(j);
	return *this;
}

MatrixWithOp::mat_mut_ptr_t
MultiVectorMutableCols::clone()
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

MatrixWithOp::mat_ptr_t
MultiVectorMutableCols::clone() const
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

void MultiVectorMutableCols::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	,const VectorWithOp& x, value_type b
	) const
{
	using AbstractLinAlgPack::dot;

	// y = b*y
	LinAlgOpPack::Vt_S(y,b);

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// y += a*M*x
		//
		// =>
		//
		// y += a * x(1) * M(:,1) + a * x(2) * M(:,2) + ...
		//
		for( size_type j = 1; j <= col_vecs_.size(); ++j )
			LinAlgOpPack::Vp_StV( y, a * x.get_ele(j), *col_vecs_[j-1] );
	}
	else {
		//
		// y += a*M'*x
		//
		// =>
		//
		// y(1) += a M(:,1)*x
		// y(2) += a M(:,2)*x
		// ...
		//
		for( size_type j = 1; j <= col_vecs_.size(); ++j )
			y->set_ele(
				j
				,y->get_ele(j) + a * dot(*col_vecs_[j-1],x)
				);
	}
}

void MultiVectorMutableCols::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	,const SpVectorSlice& x, value_type b
	) const
{
	using AbstractLinAlgPack::dot;

	// y = b*y
	LinAlgOpPack::Vt_S(y,b);

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// y += a*M*x
		//
		// =>
		//
		// y += a * x(1) * M(:,1) + a * x(2) * M(:,2) + ...
		//
		SpVectorSlice::difference_type o = x.offset();
		for( SpVectorSlice::const_iterator itr = x.begin(); itr != x.end(); ++itr ) {
			const size_type j = itr->index() + o;
			LinAlgOpPack::Vp_StV( y, a * itr->value(), *col_vecs_[j-1] );
		}
	}
	else {
		//
		// y += a*M'*x
		//
		// =>
		//
		// y(1) += a M(:,1)*x
		// y(2) += a M(:,2)*x
		// ...
		//
		for( size_type j = 1; j <= col_vecs_.size(); ++j )
			y->set_ele(
				j
				,y->get_ele(j) + a * dot(*col_vecs_[j-1],x)
				);
	}
}

void MultiVectorMutableCols::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	using LinAlgOpPack::dot;
	MatrixSymWithOpGetGMSSymMutable
		*symwo_gms_lhs = dynamic_cast<MatrixSymWithOpGetGMSSymMutable*>(sym_lhs);
	if(!symwo_gms_lhs) {
		MatrixWithOp::syrk(M_trans,alpha,beta,sym_lhs); // Boot it
		return;
	}
	MatrixDenseSymMutableEncap  sym_gms(symwo_gms_lhs);
	const int num_vecs = this->col_vecs_.size();
	THROW_EXCEPTION(
		num_vecs != sym_gms().rows(), std::logic_error
		,"MultiVectorMutableCols::syrk(...) : Error, sizes do not match up!" );
	// Fill the upper or lower triangular region.
	if( sym_gms().uplo() == BLAS_Cpp::upper ) {
		for( int i = 1; i <= num_vecs; ++i ) {
			for( int j = i; j <= num_vecs; ++j ) { // Upper triangle!
				sym_gms().gms()(i,j) = dot(*col_vecs_[i-1],*col_vecs_[j-1]);
			}
		}
	}
	else {
		for( int i = 1; i <= num_vecs; ++i ) {
			for( int j = 1; j <= i; ++j ) { // Lower triangle!
				sym_gms().gms()(i,j) = dot(*col_vecs_[i-1],*col_vecs_[j-1]);
			}
		}
	}
}

// Overridden from MultiVector

MultiVector::access_by_t
MultiVectorMutableCols::access_by() const
{
	return COL_ACCESS;
}

// Overridden from MultiVectorMutable

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableCols::col(index_type j)
{
	THROW_EXCEPTION( !(  1 <= j  && j <= col_vecs_.size() ), std::logic_error, "Error!" );
	return col_vecs_[j-1];
}

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableCols::row(index_type i)
{
	return MemMngPack::null;
}

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableCols::diag(int k)
{
	return MemMngPack::null;
}

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutableCols::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
#ifdef _DEBUG
	const size_type rows = this->rows();
	THROW_EXCEPTION(
		!( row_rng.full_range() || (row_rng.lbound() == 1 && row_rng.ubound() == rows) )
		,std::logic_error, "Error, can't handle subrange on vectors yet!" );
#endif
	return MemMngPack::rcp(
		new MultiVectorMutableCols(
			space_cols_,space_rows_->sub_space(col_rng),&col_vecs_[col_rng.lbound()-1]
			) );
}
	
} // end namespace SparseLinAlgPack
