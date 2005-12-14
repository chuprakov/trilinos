// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>

#include <algorithm>

#include "AbstractLinAlgPack_ConvertToCSC.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Mp_StM;
}

namespace AbstractLinAlgPack {

void ConvertToSparseCompressedColumnPack::vector_insert_nonzeros(
	  const DVectorSlice&				vs
	, value_type						alpha
	, size_type							row_offset
	, size_type							col_j
	, const IVector::value_type*		row_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	size_type ele = next_nz_in_col[ col_j - 1 ];
	next_nz_in_col[ col_j - 1 ] += vs.size();
	// Input the nonzero elements in D
	BLAS_Cpp::copy( vs.size(), vs.start_ptr(), vs.stride(), D_val + (ele - 1), 1 );
	BLAS_Cpp::scal( vs.size(), alpha, D_val + (ele - 1), 1 );
	// Set the row indices
	if(D_row_i) {
		D_row_i += (ele - 1);
		for(size_type i = 1; i <= vs.size(); ++i )
			*D_row_i++ = row_perm[ row_offset + i - 1 ];		
	}
}

size_type ConvertToSparseCompressedColumnPack::dense_num_in_column(
	  size_type							rows
	, size_type							cols
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	)
{
	if(trans == BLAS_Cpp::trans)
		std::swap(rows,cols);
	for( size_type j = 1; j <= cols; ++j )
		num_in_col[ col_perm[ col_offset + j - 1 ] - 1 ] += rows;
	return rows * cols;
}

void ConvertToSparseCompressedColumnPack::dense_insert_nonzeros(
	  const DMatrixSlice&				gms
	, BLAS_Cpp::Transp					trans
	, value_type						alpha
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	using DenseLinAlgPack::col;
	size_type cols = ( trans == BLAS_Cpp::no_trans ? gms.cols() : gms.rows() );
	for( size_type j = 1; j <= cols; ++j) {
		vector_insert_nonzeros( col( gms, trans, j ), alpha, row_offset
			, col_perm[ col_offset + j - 1 ], row_perm, next_nz_in_col, D_val, D_row_i );
	}
}

value_type ConvertToSparseCompressedColumnPack::dense_insert_scaled_nonzeros(
	  const DMatrixSlice&				gms
	, BLAS_Cpp::Transp					trans
	, value_type						scaled_max_ele
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	using DenseLinAlgPack::norm_inf;
	value_type alpha = 0;
	for( size_type j = 1; j <= gms.cols(); ++j )
		alpha = std::_MAX( alpha, norm_inf( gms.col(j) ) );
	// scaled_max_ele = max|alpha*A| = alpha * max|A| 
	alpha = scaled_max_ele / alpha;
	dense_insert_nonzeros( gms, trans, alpha, row_offset, col_offset, row_perm
		, col_perm, next_nz_in_col, D_val, D_row_i );
	return alpha;
}

size_type ConvertToSparseCompressedColumnPack::num_in_column(
	  const MatrixOp&				m
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	)
{
	const ConvertToCSC*
		conv_m = dynamic_cast<const ConvertToCSC*>( &m );
	if(conv_m)
		return conv_m->num_in_column( trans, col_offset, col_perm, num_in_col );
	else
		return dense_num_in_column( m.rows(), m.cols(), trans, col_offset, col_perm, num_in_col );		 
}

void ConvertToSparseCompressedColumnPack::insert_nonzeros(
	  const MatrixOp&				m
	, BLAS_Cpp::Transp					trans
	, value_type						alpha
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	using LinAlgOpPack::assign;
	const ConvertToCSC*
		conv_m = dynamic_cast<const ConvertToCSC*>( &m );
	if(conv_m) {
		conv_m->insert_nonzeros( trans, alpha, row_offset, col_offset, row_perm
			, col_perm, next_nz_in_col, D_val, D_row_i );
	}
	else {
		DMatrix _m;
		assign( &_m, m, BLAS_Cpp::no_trans );
		dense_insert_nonzeros( _m, trans, alpha, row_offset, col_offset, row_perm
			, col_perm, next_nz_in_col, D_val, D_row_i );		 
	}
}

value_type ConvertToSparseCompressedColumnPack::insert_scaled_nonzeros(
	  const MatrixOp&				m
	, BLAS_Cpp::Transp					trans
	, value_type						scaled_max_ele
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	using LinAlgOpPack::assign;
	const ConvertToCSC*
		conv_m = dynamic_cast<const ConvertToCSC*>( &m );
	if(conv_m) {
		return conv_m->insert_scaled_nonzeros( trans, scaled_max_ele, row_offset
			, col_offset, row_perm, col_perm, next_nz_in_col, D_val, D_row_i );
	}
	else {
		DMatrix _m;
		assign( &_m, m, BLAS_Cpp::no_trans );
		return dense_insert_scaled_nonzeros( _m, trans, scaled_max_ele, row_offset
			, col_offset, row_perm, col_perm, next_nz_in_col, D_val, D_row_i );		 
	}
}

}	// end namespace AbstractLinAlgPack
