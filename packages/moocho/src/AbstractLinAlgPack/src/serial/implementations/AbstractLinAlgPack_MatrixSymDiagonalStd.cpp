// ////////////////////////////////////////////////////////////////
// MatrixSymDiagonalStd.cpp
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

#include "SparseLinAlgPack/src/MatrixSymDiagonalStd.h"
#include "SparseLinAlgPack/src/SpVectorClass.h"
#include "LinAlgPack/src/VectorOp.h"
#include "LinAlgPack/src/GenMatrixOp.h"
#include "LinAlgPack/src/LinAlgPackAssertOp.h"

namespace SparseLinAlgPack {

MatrixSymDiagonalStd::MatrixSymDiagonalStd()
{}

VectorSlice MatrixSymDiagonalStd::diag()
{
	return diag_();
}

const VectorSlice MatrixSymDiagonalStd::diag() const
{
	return diag_();
}

// Overridden from MatrixSymInitDiagonal

void MatrixSymDiagonalStd::init_identity( size_type n, value_type alpha = 1.0 )
{
	diag_.resize(n);
	if(n)
		diag_ = alpha;
}

void MatrixSymDiagonalStd::init_diagonal( const VectorSlice& diag )
{
	diag_ = diag;
}

// Overridden from Matrix

size_type MatrixSymDiagonalStd::rows() const
{
	return diag_.size();
}

// Overridden from MatrixWithOp

void MatrixSymDiagonalStd::Mp_StM(
	GenMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
	using LinAlgPack::Vp_StV;
	// Assert that the dimensions match up.
	LinAlgPack::Mp_M_assert_sizes( gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), trans_rhs );
	// Add to the diagonal
	Vp_StV( &gms_lhs->diag(), alpha, diag_ );
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2, value_type beta) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	//
	// y = b*y + a * op(A) * x
	//
	LinAlgPack::Vp_MtV_assert_sizes(
		vs_lhs->size(), n, n, trans_rhs1, vs_rhs2.size() );
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) += a * diag(j) * x(j), for j = 1...n
	//
	if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
		// Optimized implementation
		const value_type
			*d_itr      = diag.raw_ptr(),
			*x_itr      = vs_rhs2.raw_ptr();
		value_type
			*y_itr      = vs_lhs->raw_ptr(),
			*y_end      = y_itr + vs_lhs->size();

		if( beta == 0.0 ) {
			while( y_itr != y_end )
				*y_itr++ = alpha * (*d_itr++) * (*x_itr++);
		}
		else if( beta == 1.0 ) {
			while( y_itr != y_end )
				*y_itr++ += alpha * (*d_itr++) * (*x_itr++);
		}
		else {
			for( ; y_itr != y_end; ++y_itr )
				*y_itr = beta * (*y_itr) + alpha * (*d_itr++) * (*x_itr++);
		}
	}
	else {
		// Generic implementation
		VectorSlice::const_iterator
			d_itr = diag.begin(),
			x_itr = vs_rhs2.begin();
		VectorSlice::iterator
			y_itr = vs_lhs->begin(),
			y_end = vs_lhs->end();
		for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
#ifdef LINALGPACK_CHECK_RANGE
			assert( d_itr < diag.end() );
			assert( x_itr < vs_rhs2.end() );
			assert( y_itr < vs_lhs->end() );
#endif
			*y_itr = beta * (*y_itr) + alpha * (*d_itr) * (*x_itr);
		}
	}
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, sv_rhs2.size() );
	//
	// y = b*y + a * op(A) * x
	//
	LinAlgPack::Vt_S(vs_lhs,beta); // y = b * y
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) += a * diag(j) * x(j), for j = 1...n
	//
	// x is sparse so take account of this.

	for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
		 ; x_itr != sv_rhs2.end()
		 ; ++x_itr )
	{
		(*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
			+= alpha * diag(x_itr->indice() + sv_rhs2.offset()) * x_itr->value();
			// Note: The indice x(i) invocations are ranged check
			// if this is compiled into the code.
	}
}

// Overridden from MatrixWithOpFactorized

void MatrixSymDiagonalStd::V_InvMtV(
	VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = inv(op(A)) * x
	//
	// A is symmetric and diagonal (A = diag(diag)) so:
	//
	// y(j) = x(j) / diag(j), for j = 1...n

	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, vs_rhs2.size() );
	
	if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
		// Optimized implementation
		const value_type
			*d_itr      = diag.raw_ptr(),
			*x_itr      = vs_rhs2.raw_ptr();
		value_type
			*y_itr      = vs_lhs->raw_ptr(),
			*y_end      = y_itr + vs_lhs->size();
		while( y_itr != y_end )
			*y_itr++ = (*x_itr++) / (*d_itr++);
	}
	else {
		// Generic implementation
		VectorSlice::const_iterator
			d_itr = diag.begin(),
			x_itr = vs_rhs2.begin();
		VectorSlice::iterator
			y_itr = vs_lhs->begin(),
			y_end = vs_lhs->end();
		for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
			assert( d_itr < diag.end() );
			assert( x_itr < vs_rhs2.end() );
			assert( y_itr < vs_lhs->end() );
			*y_itr = (*x_itr)/(*d_itr);
		}
	}
}

void MatrixSymDiagonalStd::V_InvMtV(
	VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = inv(op(A)) * x
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) = x(j) / diag(j), for j = 1...n
	//
	// x is sparse so take account of this.
	
	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, sv_rhs2.size() );

	for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
		 ; x_itr != sv_rhs2.end()
		 ; ++x_itr )
	{
		(*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
			= x_itr->value() / diag(x_itr->indice() + sv_rhs2.offset());
			// Note: The indice x(i) invocations are ranged check
			// if this is compiled into the code.
	}
}

} // end namespace SparseLinAlgPack
