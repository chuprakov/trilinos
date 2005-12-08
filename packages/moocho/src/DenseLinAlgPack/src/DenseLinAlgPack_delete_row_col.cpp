// ///////////////////////////////////////////////////////
// delete_row_col.cpp
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

#include "DenseLinAlgPack_delete_row_col.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

void DenseLinAlgPack::delete_row_col( size_type kd, DMatrixSliceTriEle* tri_M )
{
	// Validate input
	assert( tri_M );
	assert( tri_M->rows() );
	assert( 1 <= kd && kd <= tri_M->rows() );

	DMatrixSlice   M = tri_M->gms();
	const size_type  n = M.rows();

	if( tri_M->uplo() == BLAS_Cpp::lower ) {
		// Move M31 up one row at a time
		if( 1 < kd && kd < n ) {
			Range1D rng(1,kd-1);
			for( size_type i = kd; i < n; ++i )
				M.row(i)(rng) = M.row(i+1)(rng);
		}
		// Move M33 up and to the left one column at a time
		if( kd < n ) {
			for( size_type i = kd; i < n; ++i )
				M.col(i)(i,n-1) = M.col(i+1)(i+1,n);
		}
	}
	else if(  tri_M->uplo() == BLAS_Cpp::upper ) {
		// Move M13 left one column at a time.
		if( 1 < kd && kd < n ) {
			Range1D rng(1,kd-1);
			for( size_type j = kd; j < n; ++j )
				M.col(j)(rng) = M.col(j+1)(rng);
		}
		// Move the updated U33 up and left one column at a time.
		if( kd < n ) {
			for( size_type j = kd; j < n; ++j )
				M.col(j)(kd,j) = M.col(j+1)(kd+1,j+1);
		}
	}
	else {
		assert(0); // Invalid input
	}
}
