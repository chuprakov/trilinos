// /////////////////////////////////////////////////////////////////
// MatrixSymDiagonalSparse.cpp
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

#include <fstream>		// For debugging only

#include "SparseLinAlgPack/include/MatrixSymDiagonalSparse.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/EtaVector.h"
#include "AbstractLinAlgPack/include/AbstractLinAlgPackAssertOp.h"
#include "SparseLinAlgPack/include/SpVectorOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/GenMatrixAssign.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/GenMatrixOp.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace SparseLinAlgPack {

MatrixSymDiagonalSparse::MatrixSymDiagonalSparse()
	: num_updates_at_once_(0)	// Flag that it is to be determined internally.
{}

// Overridden from MatrixBase

size_type MatrixSymDiagonalSparse::rows() const
{
	return diag().dim();
}

// Overridden from MatrixWithOp

std::ostream& MatrixSymDiagonalSparse::output(std::ostream& out) const
{
	out	<< "*** Sparse diagonal matrix ***\n"
		<< "diag =\n" << diag();
	return out;
}

// Overridden from MatrixWithOpSerial

void MatrixSymDiagonalSparse::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta) const
{
	const SpVectorSlice &diag = this->diag();

	size_type n = diag.dim();

	// Assert that the dimensions of the aruments match up and if not
	// then throw an excption.
	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->dim(), n, n, trans_rhs1, vs_rhs2.dim() );

	// y = b*y + a * op(A) * x
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) = b*y(j) + a * diag(j) * x(j), for j = 1...n

	for( SpVectorSlice::const_iterator d_itr = diag.begin(); d_itr != diag.end(); ++d_itr )
	{
		const size_type i = d_itr->index(); 
		(*vs_lhs)(i) = beta * (*vs_lhs)(i) + alpha * d_itr->value() * vs_rhs2(i);
	}
}

// Overridden from MatrixSymWithOpSerial

void MatrixSymDiagonalSparse::Mp_StMtMtM(
	sym_gms* B, value_type alpha
	,EMatRhsPlaceHolder dummy_place_holder
	,const MatrixWithOpSerial& A, BLAS_Cpp::Transp A_trans
	,value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using BLAS_Cpp::trans_not;

	using LinAlgPack::nonconst_tri_ele;
	using LinAlgPack::assign;
	using LinAlgPack::syrk;
	using LinAlgPack::assert_print_nan_inf;

	using LinAlgOpPack::V_MtV;

	typedef EtaVector eta;

	// Assert the size matching of M * op(A)
	LinAlgPack::MtV_assert_sizes(
		  this->rows(), this->cols(), BLAS_Cpp::no_trans
		, rows( A.rows(), A.cols(), A_trans ) );

	// Assert size matchin of B = (op(A') * M * op(A))
	LinAlgPack::Vp_V_assert_sizes(
		  B->cols(), cols( A.rows(), A.cols(), A_trans ) );

	//
	// B = a * op(A') * M * op(A)
	//
	//   = a * op(A') * M^(1/2) * M^(1/2) * op(A)
	//
	//   = a * E * E'
	//
	// E = M^(1/2) * op(A)
	//
	//     [ .                                                 ] [ .              ]
	//     [   sqrt(M(j(1)))                                   ] [ op(A)(j(1),:)  ]
	//     [                .                                  ] [ .              ]
	//   = [                  sqrt(M(j(i))                     ] [ op(A)(j(i),:)  ]
	//     [                              .                    ] [ .              ]
	//     [                                sqrt(M(j(nz))      ] [ op(A)(j(nz),:) ]
	//     [                                               .   ] [ .              ]
	//
	//
	//     [ .        ]
	//     [ d(j(1))' ]
	//     [ .        ]
	//   = [ d(j(i))' ]
	//     [ .        ]
	//     [ d(j(1))' ]
	//     [ .        ]
	//
	//     where: d(j(i)) = sqrt(M(j(i)) * op(A')(:,j(i))    <: R^m 
	//                    = sqrt(M(j(i)) * op(A') * e(j(i))  <: R^m
	//
	//  Above M^(1/2) only has nz nonzero elements sqrt(M(j(i)), i = 1..nz and only
	//  the corresponding rows of op(A)(j(i),:), i = 1..nz are shown.  A may in fact
	//  dense matrix but the columns are extracted through op(A)*eta(j(i)), i=1..nz.
	//
	//  The above product B = a * E * E' is a set of nz rank-1 updates and can be written
	//  in the form:
	//
	//  B = sum( a * d(j(i)) * d(j(i))', i = 1..nz )
	//
	//  Since it is more efficient to perform several rank-1 updates at a time we will
	//  perform them in blocks.
	//
	//  B = B + D(k) * D(k)', k = 1 .. num_blocks
	//
	//      where:
	//         num_blocks = nz / num_updates_at_once + 1 (integer division)
	//         D(k) = [ d(j(i1)) ... d(j(i2)) ]
	//         i1 = (k-1) * num_updates_at_once + 1
	//         i2 = i1 + num_updates_at_once - 1
	
	const SpVectorSlice
		&diag = this->diag();

	const size_type
		n = this->rows(),
		m = cols(A.rows(),A.cols(),A_trans);

	// Get the actual number of updates to use per rank-(num_updates) update
	const size_type
		num_updates
			= std::_MIN( num_updates_at_once()
							? num_updates_at_once()
							: 20	// There may be a better default value for this?
						, diag.nz()
						);

	// Get the number of blocks of rank-(num_updates) updates
	size_type
		num_blocks = diag.nz() / num_updates;
	if( diag.nz() % num_updates > 0 )
		num_blocks++;

	// Initialize B = b*B
	if( b != 1.0 )
		assign( &nonconst_tri_ele( B->gms(), B->uplo() ), 0.0 );

	// Perform the rank-(num_updates) updates
	GenMatrix D(m,num_updates);
	for( size_type k = 1; k <= num_blocks; ++k ) {
		const size_type
			i1 = (k-1) * num_updates + 1,
			i2 = std::_MIN( diag.nz(), i1 + num_updates - 1 );
		// Generate the colunns of D(k)
		SpVectorSlice::const_iterator
			m_itr = diag.begin() + (i1-1);
		for( size_type l = 1; l <= i2-i1+1; ++l, ++m_itr ) {
			assert( m_itr < diag.end() );
			assert( m_itr->value() >= 0.0 );
			V_MtV( &D.col(l), A, trans_not(A_trans)
				, eta( m_itr->index(), n, ::sqrt(m_itr->value()) )() );
		}
		const GenMatrixSlice
			D_update = D(1,m,1,i2-i1+1);


//		// For debugging only
//		std::ofstream ofile("MatrixSymDiagonalSparse_Error.out");
// 		assert_print_nan_inf( D_update, "D", true, &ofile );
		// Perform the rank-(num_updates) update
		syrk( BLAS_Cpp::no_trans, alpha, D_update, 1.0, B );
	}
}

// Overridden from MatrixConvertToSparseFortranCompatible

FortranTypes::f_int
MatrixSymDiagonalSparse::num_nonzeros( EExtractRegion extract_region ) const
{
	return diag().nz();
}

void MatrixSymDiagonalSparse::coor_extract_nonzeros(
	EExtractRegion                extract_region
	,const FortranTypes::f_int    len_Aval
	,FortranTypes::f_dbl_prec     Aval[]
	,const FortranTypes::f_int    len_Aij
	,FortranTypes::f_int          Arow[]
	,FortranTypes::f_int          Acol[]
	,const FortranTypes::f_int    row_offset
	,const FortranTypes::f_int    col_offset
	 ) const
{
	const SpVectorSlice
		&diag = this->diag();

	assert( len_Aval == 0 || (len_Aval == diag.nz() && Aval)			);
	assert( len_Aij  == 0 || (len_Aij  == diag.nz() && Arow && Acol)	);

	if( len_Aval > 0 ) {
		SpVectorSlice::const_iterator
			itr;
		FortranTypes::f_dbl_prec
			*l_Aval;
		for( itr = diag.begin(), l_Aval = Aval; itr != diag.end(); ++itr ) {
			*l_Aval++ = itr->value();
		}			
	}

	if( len_Aij > 0 ) {
		SpVectorSlice::const_iterator
			itr;
		FortranTypes::f_int
			*l_Arow, *l_Acol;
		for( itr = diag.begin(), l_Arow = Arow, l_Acol = Acol; itr != diag.end(); ++itr ) {
			const FortranTypes::f_int
				ij = itr->index() + diag.offset();
			*l_Arow++ = ij + row_offset;
			*l_Acol++ = ij + col_offset;
		}			
	}

}

}	// end namespace SparseLinAlgPack
