// /////////////////////////////////////
// AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDef.hpp
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

#ifndef COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DEF_H
#define COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DEF_H

#include "AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDecl.hpp"

namespace AbstractLinAlgPack {

template<class T_COOM>
size_type COOM_num_in_column(
	  const T_COOM&						m
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	)
{
	if(!m.nz()) return 0;
	if( trans == BLAS_Cpp::no_trans ) {
		// non transposed
		typename T_COOM::difference_type loc_co = m.col_offset(); 
		for( typename T_COOM::const_iterator itr = m.begin(); itr != m.end(); ++itr )
			num_in_col[ col_perm[ col_offset + loc_co + itr->col_j()  -1 ]  -1 ]++;
	}
	else {
		// transposed
		typename T_COOM::difference_type loc_ro = m.row_offset(); 
		for( typename T_COOM::const_iterator itr = m.begin(); itr != m.end(); ++itr ) {
			const size_type i = itr->row_i();
			num_in_col[ col_perm[ col_offset + loc_ro + i - 1 ]  - 1 ]++;
		}
	}
	return m.nz();
}

template<class T_COOM>
void COOM_insert_nonzeros(
	  const T_COOM&						m
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
	if(!m.nz()) return;
	typename T_COOM::difference_type
		loc_ro = m.row_offset(),
		loc_co = m.col_offset();
	if( trans == BLAS_Cpp::no_trans ) {
		// non transposed
		for( typename T_COOM::const_iterator itr = m.begin(); itr != m.end(); ++itr ) {
			const size_type
				i = loc_ro + itr->row_i(),
				j = loc_co + itr->col_j();
			const size_type
				ele = next_nz_in_col[ col_perm[ col_offset + j - 1 ] - 1 ]++;
			D_val[ ele - 1 ] = alpha * itr->value();
			if(D_row_i)
				D_row_i[ ele - 1 ] = row_perm[ row_offset + i - 1 ];
		}
	}
	else {
		// transposed
		for( typename T_COOM::const_iterator itr = m.begin(); itr != m.end(); ++itr ) {
			const size_type
				i = loc_co + itr->col_j(),
				j = loc_ro + itr->row_i();
			const size_type
				ele = next_nz_in_col[ col_perm[ col_offset + j - 1 ] - 1 ]++;
			D_val[ ele - 1 ] = alpha * itr->value();
			if(D_row_i)
				D_row_i[ ele - 1 ] = row_perm[ row_offset + i - 1 ];
		}
	}
}

template<class T_COOM>
value_type COOM_insert_scaled_nonzeros(
	  const T_COOM&						m
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
	value_type alpha = 0;
	for( typename T_COOM::const_iterator itr = m.begin(); itr != m.end(); ++itr ) {
		register const value_type val = ::fabs( itr->value() );
		if( val > alpha ) alpha = val;
	}
	// scaled_max_ele = max|alpha*A| = alpha * max|A| 
	alpha = scaled_max_ele / alpha;
	COOM_insert_nonzeros( m, trans, alpha, row_offset
		, col_offset, row_perm, col_perm, next_nz_in_col, D_val, D_row_i );
	return alpha;
}

} // end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DEF_H
