// //////////////////////////////////////////////////////////////////////////
// GenMatrixOpBLAS.cpp
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

#include "../include/GenMatrixOp.h"
#include "../include/BLAS_Cpp.h"

namespace {

using LinAlgPack::Vector;
using LinAlgPack::VectorSlice;
using LinAlgPack::GenMatrix;
using LinAlgPack::GenMatrixSlice;
using LinAlgPack::col;
using LinAlgPack::value_type;
using LinAlgPack::assert_gms_sizes;
using LinAlgPack::tri_ele_gms;
using LinAlgPack::tri_gms;
using LinAlgPack::sym_gms;
using LinAlgPack::assign;
using BLAS_Cpp::Transp;
using BLAS_Cpp::no_trans;
using BLAS_Cpp::trans;
using BLAS_Cpp::Uplo;
using BLAS_Cpp::upper;
using BLAS_Cpp::lower;

}	// end namespace

// ///////////////////////////////////////////////////////////////////////////////////
// Assignment Fucntions

namespace {

// implementation: gms_lhs = alpha (elementwise)
inline void i_assign(GenMatrixSlice* gms_lhs, value_type alpha) {
	for(GenMatrixSlice::size_type i = 1; i <= gms_lhs->cols(); ++i)
		gms_lhs->col(i) = alpha;	
}

// implementaion: gms_lhs = gms_rhs. Most basic copy function for rectangular matrices.
inline void i_assign_basic(GenMatrixSlice* gms_lhs, const GenMatrixSlice& gms_rhs
	,  BLAS_Cpp::Transp trans_rhs)
{
	for(GenMatrixSlice::size_type i = 1; i <= gms_lhs->cols(); ++i)
			gms_lhs->col(i) = col(gms_rhs,trans_rhs,i);
}

// implementaion: gms_lhs = op(gms_rhs). Checks for overlap and creates temporaries accordingly.
inline void i_assign(GenMatrixSlice* gms_lhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs)
{
	switch(gms_lhs->overlap(gms_rhs)) {
		case LinAlgPack::NO_OVERLAP: // no overlap so just perform the copy
			i_assign_basic(gms_lhs,gms_rhs,trans_rhs);
			return;
		case LinAlgPack::SAME_MEM:
			if(trans_rhs == BLAS_Cpp::no_trans) return; // assignment to self, nothing to do.
		default: // either same memory that needs to be transposed or some overlap so just generate temp.
			GenMatrix temp = gms_rhs;
			i_assign_basic(gms_lhs,temp,trans_rhs);
			return;
	}
}

// Copy one triangular region into another. Does not check sizes or aliasing of argument matrices.
// A row of a upper triangular region corresponds to a col of a BLAS_Cpp::lower triangular region.
inline void i_assign_basic(tri_ele_gms* tri_lhs, const tri_ele_gms& tri_rhs)
{
	GenMatrixSlice::size_type n = tri_lhs->gms().cols();

	// Access BLAS_Cpp::lower tri by col and upper tri by row
	BLAS_Cpp::Transp
		trans_lhs = (tri_lhs->uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans,
		trans_rhs = (tri_rhs.uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans;
	
	for(int i = 1; i <= n; ++i) { // Only copy the part of the vec in tri region.
		col(tri_lhs->gms(),trans_lhs,i)(i,n) = col(tri_rhs.gms(),trans_rhs,i)(i,n);
	}
}

}	// end namespace

// gm_lhs = alpha (elementwise)
void LinAlgPack::assign(GenMatrix* gm_lhs, value_type alpha)
{
	if(!gm_lhs->rows()) gm_lhs->resize(1,1,alpha);
	i_assign(&(*gm_lhs)(),alpha);
}

// gm_lhs = op(gms_rhs)
void LinAlgPack::assign(GenMatrix* gm_lhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs)
{
	if(gm_lhs->overlap(gms_rhs) == SAME_MEM && trans_rhs == BLAS_Cpp::no_trans) return;	// assignment to self
	if(gm_lhs->overlap(gms_rhs) != NO_OVERLAP) {
		// some overlap so we must create a copy
		GenMatrix tmp(gms_rhs);
		resize_gm_lhs(gm_lhs,gms_rhs.rows(),gms_rhs.cols(),trans_rhs);
		i_assign(&(*gm_lhs)(), tmp, trans_rhs);
	}
	else {
		// no overlap so just assign
		resize_gm_lhs(gm_lhs,gms_rhs.rows(),gms_rhs.cols(),trans_rhs);
		i_assign(&(*gm_lhs)(), gms_rhs, trans_rhs);
	}
}

// gms_lhs = alpha (elementwise)
void LinAlgPack::assign(GenMatrixSlice* gms_lhs, value_type alpha)
{
	if(!gms_lhs->rows())
		throw std::length_error("assign(gms_lhs,alpha):  You can not assign"
								" a scalar to an unsized GenMatrixSlice" );
	else i_assign(gms_lhs, alpha);
}

// gms_lhs = op(gms_rhs)
void LinAlgPack::assign(GenMatrixSlice* gms_lhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs)
{
	assert_gms_lhs(*gms_lhs,gms_rhs.rows(),gms_rhs.cols(),trans_rhs);
	i_assign(gms_lhs, gms_rhs, trans_rhs);
}

// tri_lhs = alpha (elementwise)
void LinAlgPack::assign(tri_ele_gms* tri_lhs, value_type alpha)
{
	if(!tri_lhs->gms().rows())
		throw std::length_error("assign(tri_lhs,alpha):  You can not assign"
								" a scalar to an unsized triangular GenMatrixSlice" );
	assert_gms_square(tri_lhs->gms()); 
	GenMatrixSlice::size_type n = tri_lhs->gms().cols();
	// access BLAS_Cpp::lower tri by col and upper tri by row
	BLAS_Cpp::Transp
		trans_lhs = (tri_lhs->uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans;
	for(int i = 1; i <= n; ++i)
		col( tri_lhs->gms(), trans_lhs , i )(i,n) = alpha;
}

// tri_lhs = tri_rhs
void LinAlgPack::assign(tri_ele_gms* tri_lhs, const tri_ele_gms& tri_rhs)
{
	assert_gms_lhs(tri_lhs->gms(),tri_rhs.gms().rows(),tri_rhs.gms().cols(),BLAS_Cpp::no_trans);
		
	switch(tri_lhs->gms().overlap(tri_rhs.gms())) {
		case SAME_MEM:
			if(tri_lhs->uplo() == tri_rhs.uplo()) return; // Assignment to self to exit

		case SOME_OVERLAP:
			// Test for the special case where the upper tri region (above diagonal) of a
			// matrix is being copied into the BLAS_Cpp::lower tri region (below the diagonal) of
			// the same matrix or visa-versa.  No temp is need in this case
			if(tri_lhs->uplo() != tri_rhs.uplo()) {
				const GenMatrixSlice	*up = (tri_lhs->uplo() == upper)
												? &tri_lhs->gms() : &tri_rhs.gms();
				const GenMatrixSlice	*lo = (tri_rhs.uplo() == BLAS_Cpp::lower)
												? &tri_lhs->gms() : &tri_rhs.gms();
				if(lo->col_ptr(1) + lo->max_rows() - 1 == up->col_ptr(1)) { // false if SAME_MEM
					// One triangular region is being copied into another so no temp needed.
					i_assign_basic(tri_lhs, tri_rhs);
					return; 
				}
			}
			// Give up and copy the vs_rhs as a temp.
			{
				GenMatrix temp(tri_rhs.gms());
				i_assign_basic(tri_lhs, tri_ele(temp,tri_rhs.uplo()));
				return;
			}

		case NO_OVERLAP:	// no overlap so just perform the copy
			i_assign_basic(tri_lhs, tri_rhs);
			return;
	}
}

// /////////////////////////////////////////////////////////////////////////////////////////
// Element-wise Algebraic Operations

// gms_lhs *= alpha (BLAS xSCAL)
void LinAlgPack::Mt_S(GenMatrixSlice* gms_lhs, value_type alpha)
{
	if(!gms_lhs->rows())
		throw std::length_error("Mt_S(gms_lhs,alpha):  You can not scale"
								" an unsized GenMatrixSlice" );
	for(int j = 1; j <= gms_lhs->cols(); ++j)
		Vt_S(&gms_lhs->col(j), alpha);
}

// tri_lhs *= alpha (BLAS xSCAL)
void LinAlgPack::Mt_S(tri_ele_gms* tri_lhs, value_type alpha) {
	if(!tri_lhs->gms().rows())
		throw std::length_error("Mt_S(tri_lhs,alpha):  You can not scale"
								" an unsized triangular GenMatrixSlice" );
	BLAS_Cpp::Transp
		trans_lhs = (tri_lhs->uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans;

	GenMatrixSlice::size_type n = tri_lhs->gms().cols();
	for(GenMatrixSlice::size_type j = 1; j <= n; ++j)
		Vt_S( &col( tri_lhs->gms(), trans_lhs , j )(j,n), alpha );	
}

// tri_ele_lhs += alpha * tri_ele_rhs (BLAS xAXPY)
void LinAlgPack::Mp_StM(tri_ele_gms* tri_lhs, value_type alpha, const tri_ele_gms& tri_ele_rhs)
{
	Mp_M_assert_sizes(tri_lhs->gms().rows(), tri_lhs->gms().cols(), BLAS_Cpp::no_trans
		, tri_ele_rhs.gms().rows(), tri_ele_rhs.gms().cols(), BLAS_Cpp::no_trans);

	BLAS_Cpp::Transp
		trans_lhs = (tri_lhs->uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans,
		trans_rhs = (tri_ele_rhs.uplo() == BLAS_Cpp::lower) ? BLAS_Cpp::no_trans : BLAS_Cpp::trans;

	GenMatrixSlice::size_type n = tri_lhs->gms().cols();
	for(GenMatrixSlice::size_type j = 1; j <= n; ++j)
		Vp_StV( &col(tri_lhs->gms(),trans_lhs,j)(j,n), alpha, col(tri_ele_rhs.gms(),trans_rhs,j)(j,n) );	
}

// LinAlgOpPack Compatable (compile time polymorphism)

// gms_lhs += alpha * op(gms_rhs) (BLAS xAXPY)
void LinAlgPack::Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, gms_rhs.rows(), gms_rhs.cols(), trans_rhs);
	for(GenMatrixSlice::size_type j = 1; j <= gms_lhs->cols(); ++j)
		Vp_StV( &gms_lhs->col(j), alpha, col(gms_rhs,trans_rhs,j) );	
}

namespace {
// Implement upper and lower copies for a symmetric matrix
// Does not check sizes.
// inline
void i_Mp_StSM( LinAlgPack::GenMatrixSlice* gms_lhs, LinAlgPack::value_type alpha
	, const LinAlgPack::tri_ele_gms& tri_ele_rhs)
{
	using LinAlgPack::Mp_StM;
	using LinAlgPack::tri_ele;
	const LinAlgPack::size_type n = gms_lhs->rows(); // same as cols
	Mp_StM( const_cast<LinAlgPack::tri_ele_gms*>( 
				&tri_ele((*gms_lhs)(2,n,1,n-1), BLAS_Cpp::lower ) )
			, alpha, tri_ele_rhs );
	Mp_StM( const_cast<LinAlgPack::tri_ele_gms*>(
				&tri_ele((*gms_lhs)(1,n-1,2,n), BLAS_Cpp::upper ) )
			, alpha, tri_ele_rhs );
}
}	// end namespace

// gms_lhs += alpha * op(sym_rhs) (BLAS xAXPY)
void LinAlgPack::Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, const sym_gms& sym_rhs
	, BLAS_Cpp::Transp trans_rhs )
{
	Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, sym_rhs.rows(), sym_rhs.cols(), trans_rhs);
	Vp_StV( &gms_lhs->diag(), alpha, sym_rhs.gms().diag() );
	const size_type n = gms_lhs->rows(); // same as cols
	switch(sym_rhs.uplo()) {
		case BLAS_Cpp::lower:
			i_Mp_StSM( gms_lhs, alpha, tri_ele(sym_rhs.gms()(2,n,1,n-1),BLAS_Cpp::lower) );
			break;
		case BLAS_Cpp::upper:
			i_Mp_StSM( gms_lhs, alpha, tri_ele(sym_rhs.gms()(1,n-1,2,n),BLAS_Cpp::upper) );
			break;
	}
}

// gms_lhs += alpha * op(tri_rhs) (BLAS xAXPY)
void LinAlgPack::Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	using BLAS_Cpp::trans;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::lower;
	using BLAS_Cpp::upper;
	Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, tri_rhs.rows(), tri_rhs.cols(), trans_rhs );
	// diagonal
	if( tri_rhs.diag() == BLAS_Cpp::nonunit )
		Vp_StV( &gms_lhs->diag(), alpha, tri_rhs.gms().diag() );
	else
		Vp_S( &gms_lhs->diag(), alpha );
	// upper or lower triangle
	const size_type n = gms_lhs->rows(); // same as cols
	if( n == 1 )
		return;	// There is not lower or upper triangular region
	const tri_ele_gms
		tri = ( tri_rhs.uplo() == lower
				?		tri_ele(tri_rhs.gms()(2,n,1,n-1),lower)
					 :  tri_ele(tri_rhs.gms()(1,n-1,2,n),upper) );
	const BLAS_Cpp::Uplo
		as_uplo = (		( tri_rhs.uplo() == lower && trans_rhs == no_trans )
					||	( tri_rhs.uplo() == upper && trans_rhs == trans )
					?	lower : upper											);
	switch(as_uplo) {
		case lower:
			Mp_StM( const_cast<LinAlgPack::tri_ele_gms*>(
						&tri_ele((*gms_lhs)(2,n,1,n-1),lower) )
					, alpha, tri );
			break;
		case upper:
			Mp_StM( const_cast<LinAlgPack::tri_ele_gms*>(
						&tri_ele((*gms_lhs)(1,n-1,2,n),upper) )
					, alpha, tri );
			break;
	}
}

// /////////////////////////////////////////////////////////////////////////////////////
// Level-2 BLAS (vector-matrtix) Liner Algebra Operations

/// vs_lhs = alpha * op(gms_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
void LinAlgPack::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta)
{
	Vp_MtV_assert_sizes(vs_lhs->size(), gms_rhs1.rows()	, gms_rhs1.cols(), trans_rhs1
		, vs_rhs2.size());
	BLAS_Cpp::gemv(trans_rhs1,gms_rhs1.rows(),gms_rhs1.cols(),alpha,gms_rhs1.col_ptr(1)
		,gms_rhs1.max_rows(), vs_rhs2.raw_ptr(),vs_rhs2.stride(),beta,vs_lhs->raw_ptr()
		,vs_lhs->stride());
}

// vs_lhs = alpha * op(sym_rhs1) * vs_rhs2 + beta * vs_lhs. (BLAS xSYMV).
void LinAlgPack::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const sym_gms& sym_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta)
{
	assert_gms_square(sym_rhs1.gms());
	Vp_MtV_assert_sizes(vs_lhs->size(), sym_rhs1.gms().rows(), sym_rhs1.gms().cols(), trans_rhs1
		, vs_rhs2.size());
	BLAS_Cpp::symv(sym_rhs1.uplo(),sym_rhs1.gms().rows(),alpha,sym_rhs1.gms().col_ptr(1)
		,sym_rhs1.gms().max_rows(),vs_rhs2.raw_ptr(),vs_rhs2.stride(),beta
		,vs_lhs->raw_ptr(),vs_lhs->stride());
}

// v_lhs = op(tri_rhs1) * vs_rhs2 (BLAS xTRMV)
void LinAlgPack::V_MtV(Vector* v_lhs, const tri_gms& tri_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	assert_gms_square(tri_rhs1.gms());
	MtV_assert_sizes(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1, vs_rhs2.size());
	v_lhs->resize(tri_rhs1.gms().rows());
	(*v_lhs) = vs_rhs2;
	BLAS_Cpp::trmv(tri_rhs1.uplo(),trans_rhs1,tri_rhs1.diag(),tri_rhs1.gms().rows()
		,tri_rhs1.gms().col_ptr(1),tri_rhs1.gms().max_rows(), v_lhs->raw_ptr(),1);
}

// vs_lhs = op(tri_rhs1) * vs_rhs2 (BLAS xTRMV)
void LinAlgPack::V_MtV(VectorSlice* vs_lhs, const tri_gms& tri_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	assert_gms_square(tri_rhs1.gms());
	MtV_assert_sizes(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1, vs_rhs2.size());
	Vp_V_assert_sizes( vs_lhs->size(), tri_rhs1.gms().rows() );
	(*vs_lhs) = vs_rhs2;
	BLAS_Cpp::trmv(tri_rhs1.uplo(),trans_rhs1,tri_rhs1.diag(),tri_rhs1.gms().rows()
		,tri_rhs1.gms().col_ptr(1),tri_rhs1.gms().max_rows(), vs_lhs->raw_ptr(),vs_lhs->stride());
}

// vs_lhs = alpha * op(tri_rhs1) * vs_rhs2 + beta * vs_lhs.
void LinAlgPack::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta)
{
	Vp_MtV_assert_sizes(vs_lhs->size(),tri_rhs1.gms().rows(),tri_rhs1.gms().cols()
						,trans_rhs1,vs_rhs2.size() );

	// If op(gms_rhs2) == gms_lhs and beta = 0.0 then this is a direct call to the BLAS.
	if( vs_lhs->overlap(vs_rhs2) == SAME_MEM && beta == 0.0 )
	{
		V_MtV(vs_lhs, tri_rhs1, trans_rhs1, vs_rhs2);
		if(alpha != 1.0) Vt_S(vs_lhs,alpha);
	}
	else {
		// This is something else so the alias problem is not handled.
		if(beta != 1.0) Vt_S(vs_lhs,beta);
		Vector tmp;
		V_MtV(&tmp, tri_rhs1, trans_rhs1, vs_rhs2);
		Vp_StV(vs_lhs,alpha,tmp());
	}
}

// v_lhs = inv(op(tri_rhs1)) * vs_rhs2 (BLAS xTRSV)
void LinAlgPack::V_InvMtV(Vector* v_lhs, const tri_gms& tri_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	assert_gms_square(tri_rhs1.gms());
	MtV_assert_sizes(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1, vs_rhs2.size());
	v_lhs->resize(tri_rhs1.gms().rows());
	(*v_lhs) = vs_rhs2;
	BLAS_Cpp::trsv(tri_rhs1.uplo(),trans_rhs1,tri_rhs1.diag(),tri_rhs1.gms().rows()
		,tri_rhs1.gms().col_ptr(1),tri_rhs1.gms().max_rows(), v_lhs->raw_ptr(),1);
}

// vs_lhs = inv(op(tri_rhs1)) * vs_rhs2 (BLAS xTRSV)
void LinAlgPack::V_InvMtV(VectorSlice* vs_lhs, const tri_gms& tri_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2)
{
	assert_gms_square(tri_rhs1.gms());
	MtV_assert_sizes(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1, vs_rhs2.size());
	Vp_V_assert_sizes( vs_lhs->size(), tri_rhs1.gms().rows() );
	(*vs_lhs) = vs_rhs2;
	BLAS_Cpp::trsv(tri_rhs1.uplo(),trans_rhs1,tri_rhs1.diag(),tri_rhs1.gms().rows()
		,tri_rhs1.gms().col_ptr(1),tri_rhs1.gms().max_rows(), vs_lhs->raw_ptr(),vs_lhs->stride());
}


// gms_lhs = alpha * vs_rhs1 * vs_rhs2' + gms_lhs (BLAS xGER)
void LinAlgPack::ger(
	value_type alpha, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2
	, GenMatrixSlice* gms_lhs )
{
	Vp_MtV_assert_sizes( vs_rhs2.size(),  gms_lhs->rows(), gms_lhs->cols()
		, BLAS_Cpp::no_trans, vs_rhs1.size() );
	BLAS_Cpp::ger(
		gms_lhs->rows(), gms_lhs->cols(), alpha
		,vs_rhs1.raw_ptr(), vs_rhs1.stride()
		,vs_rhs2.raw_ptr(), vs_rhs2.stride()
		,gms_lhs->col_ptr(1), gms_lhs->max_rows() );
}

// sym_lhs = alpha * vs_rhs * vs_rhs' + sym_lhs (BLAS xSYR).
void LinAlgPack::syr(value_type alpha, const VectorSlice& vs_rhs, sym_gms* sym_lhs)
{
	assert_gms_square(sym_lhs->gms());
	MtV_assert_sizes( sym_lhs->gms().rows(), sym_lhs->gms().cols()
		, BLAS_Cpp::no_trans, vs_rhs.size() );
	BLAS_Cpp::syr( sym_lhs->uplo(), vs_rhs.size(), alpha, vs_rhs.raw_ptr()
		, vs_rhs.stride(), sym_lhs->gms().col_ptr(1), sym_lhs->gms().max_rows() );
}

// sym_lhs = alpha * vs_rhs1 * vs_rhs2' + alpha * vs_rhs2 * vs_rhs1' + sym_lhs (BLAS xSYR2).
void LinAlgPack::syr2(value_type alpha, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2
	, sym_gms* sym_lhs)
{
	assert_gms_square(sym_lhs->gms());
	VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	MtV_assert_sizes( sym_lhs->gms().rows(), sym_lhs->gms().cols()
		, BLAS_Cpp::no_trans, vs_rhs1.size() );
	BLAS_Cpp::syr2( sym_lhs->uplo(), vs_rhs1.size(), alpha, vs_rhs1.raw_ptr()
		, vs_rhs1.stride(), vs_rhs2.raw_ptr(), vs_rhs2.stride()
		, sym_lhs->gms().col_ptr(1), sym_lhs->gms().max_rows() );
}

// //////////////////////////////////////////////////////////////////////////////////////////
// Level-3 BLAS (matrix-matrix) Linear Algebra Operations

// ////////////////////////////
// Rectangular Matrices

// gms_lhs = alpha * op(gms_rhs1) * op(gms_rhs2) + beta * gms_lhs (BLAS xGEMV).
void LinAlgPack::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
							, gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2);
	BLAS_Cpp::gemm(trans_rhs1,trans_rhs2,gms_lhs->rows(),gms_lhs->cols()
		,cols(gms_rhs1.rows(),gms_rhs1.cols(),trans_rhs1)
		,alpha,gms_rhs1.col_ptr(1),gms_rhs1.max_rows()
		,gms_rhs2.col_ptr(1),gms_rhs2.max_rows()
		,beta,gms_lhs->col_ptr(1),gms_lhs->max_rows() );
}

// ////////////////////////////
// Symmetric Matrices

namespace {

// implementation:	gms_lhs = alpha * sym_rhs * gms_rhs + beta * gms_lhs (left) (BLAS xSYMM).
// or				gms_lhs = alpha * gms_rhs * sym_rhs + beta * gms_lhs (right).
// does not check sizes.
void i_symm(BLAS_Cpp::Side side, value_type alpha, const sym_gms& sym_rhs
	, const GenMatrixSlice& gms_rhs, value_type beta, GenMatrixSlice* gms_lhs)
{
	BLAS_Cpp::symm(side,sym_rhs.uplo(),gms_lhs->rows(),gms_lhs->cols()
		,alpha,sym_rhs.gms().col_ptr(1),sym_rhs.gms().max_rows()
		,gms_rhs.col_ptr(1),gms_rhs.max_rows()
		,beta,gms_lhs->col_ptr(1),gms_lhs->max_rows() );
}

}	// end namespace

// gms_lhs = alpha * op(sym_rhs1) * op(gms_rhs2) + beta * gms_lhs (left) (BLAS xSYMM).
void LinAlgPack::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const sym_gms& sym_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, sym_rhs1.gms().rows(), sym_rhs1.gms().cols(), trans_rhs1
							, gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2);
	if(trans_rhs2 == BLAS_Cpp::no_trans) {
		i_symm(BLAS_Cpp::left,alpha,sym_rhs1,gms_rhs2,beta,gms_lhs);
	}
	else {
		// must make a temporary copy to call the BLAS
		GenMatrix tmp;
		assign(&tmp,gms_rhs2,trans_rhs2);
		i_symm(BLAS_Cpp::left,alpha,sym_rhs1,tmp(),beta,gms_lhs);
	}
}

// gms_lhs = alpha * op(gms_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (BLAS xSYMM).
void LinAlgPack::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const sym_gms& sym_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
							, sym_rhs2.gms().rows(), sym_rhs2.gms().cols(), trans_rhs2 );
	if(trans_rhs1 == BLAS_Cpp::no_trans) {
		i_symm(BLAS_Cpp::right,alpha,sym_rhs2,gms_rhs1,beta,gms_lhs);
	}
	else {
		// must make a temporary copy to call the BLAS
		GenMatrix tmp;
		assign(&tmp,gms_rhs1,trans_rhs1);
		i_symm(BLAS_Cpp::right,alpha,sym_rhs2,tmp(),beta,gms_lhs);
	}
}

// sym_lhs = alpha * op(gms_rhs) * op(gms_rhs')  + beta * sym_lhs (BLAS xSYRK).
void LinAlgPack::syrk(BLAS_Cpp::Transp trans, value_type alpha, const GenMatrixSlice& gms_rhs
	, value_type beta, sym_gms* sym_lhs)
{
	Mp_MtM_assert_sizes(	  sym_lhs->gms().rows(), sym_lhs->gms().cols(), BLAS_Cpp::no_trans
							, gms_rhs.rows(), gms_rhs.cols(), trans
							, gms_rhs.rows(), gms_rhs.cols(), trans_not(trans) );
	BLAS_Cpp::syrk(sym_lhs->uplo(),trans,sym_lhs->gms().rows()
		,cols(gms_rhs.rows(), gms_rhs.cols(), trans),alpha
		,gms_rhs.col_ptr(1),gms_rhs.max_rows(),beta
		,sym_lhs->gms().col_ptr(1),sym_lhs->gms().max_rows() );
}

// sym_lhs = alpha * op(gms_rhs1) * op(gms_rhs2') + alpha * op(gms_rhs2) * op(gms_rhs1')
void LinAlgPack::syr2k(BLAS_Cpp::Transp trans,value_type alpha, const GenMatrixSlice& gms_rhs1
	, const GenMatrixSlice& gms_rhs2, value_type beta, sym_gms* sym_lhs)
{
	Mp_MtM_assert_sizes(	  sym_lhs->gms().rows(), sym_lhs->gms().cols(), BLAS_Cpp::no_trans
							, gms_rhs1.rows(), gms_rhs1.cols(), trans
							, gms_rhs2.rows(), gms_rhs2.cols(), trans_not(trans) );
	BLAS_Cpp::syr2k(sym_lhs->uplo(),trans,sym_lhs->gms().rows()
		,cols(gms_rhs1.rows(), gms_rhs1.cols(), trans),alpha
		,gms_rhs1.col_ptr(1),gms_rhs1.max_rows()
		,gms_rhs2.col_ptr(1),gms_rhs2.max_rows(),beta
		,sym_lhs->gms().col_ptr(1),sym_lhs->gms().max_rows() );
}

// ////////////////////////////
// Triangular Matrices

// ToDo: Finish the definitions below.

namespace {

// implementation:	gms_lhs = alpha * op(tri_rhs) * gms_lhs (left) (BLAS xTRMM).
// or				gms_lhs = alpha * gms_rhs * op(tri_rhs) (right).
// does not check sizes.
inline
void i_trmm(BLAS_Cpp::Side side, BLAS_Cpp::Transp trans, value_type alpha, const tri_gms& tri_rhs
	, GenMatrixSlice* gms_lhs)
{
	BLAS_Cpp::trmm(side,tri_rhs.uplo(),trans,tri_rhs.diag(),gms_lhs->rows(),gms_lhs->cols()
		,alpha,tri_rhs.gms().col_ptr(1),tri_rhs.gms().max_rows()
		,gms_lhs->col_ptr(1),gms_lhs->max_rows() );
}

// implementation:	gms_lhs = alpha * op(tri_rhs) * op(gms_rhs) (left)
// or				gms_lhs = alpha * op(gms_rhs) * op(tri_rhs) (right)
// Takes care of temporaries but does not check sizes.
inline
void i_trmm_alt( BLAS_Cpp::Side side, value_type alpha, const tri_gms& tri_rhs
	, BLAS_Cpp::Transp trans_tri_rhs, const GenMatrixSlice& gms_rhs
	, BLAS_Cpp::Transp trans_gms_rhs, GenMatrixSlice* gms_lhs )
{
	assign(gms_lhs,gms_rhs,trans_gms_rhs);
	i_trmm( side, trans_tri_rhs, alpha, tri_rhs, gms_lhs );
}

// implementation:	gms_lhs = alpha * inv(op(tri_rhs)) * gms_lhs (left) (BLAS xTRSM).
// or				gms_lhs = alpha * gms_rhs * inv(op(tri_rhs)) (right).
// does not check sizes.
inline
void i_trsm(BLAS_Cpp::Side side, BLAS_Cpp::Transp trans, value_type alpha, const tri_gms& tri_rhs
	, GenMatrixSlice* gms_lhs)
{
	BLAS_Cpp::trsm(side,tri_rhs.uplo(),trans,tri_rhs.diag(),gms_lhs->rows(),gms_lhs->cols()
		,alpha,tri_rhs.gms().col_ptr(1),tri_rhs.gms().max_rows()
		,gms_lhs->col_ptr(1),gms_lhs->max_rows() );
}

// implementation:	gms_lhs = alpha * inv(op(tri_rhs)) * op(gms_rhs) (left)
// or				gms_lhs = alpha * op(gms_rhs) * inv(op(tri_rhs)) (right)
// Takes care of temporaries but does not check sizes.
inline
void i_trsm_alt(BLAS_Cpp::Side side, value_type alpha, const tri_gms& tri_rhs
	, BLAS_Cpp::Transp trans_tri_rhs, const GenMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_gms_rhs
	, GenMatrixSlice* gms_lhs)
{
	assign(gms_lhs,gms_rhs,trans_gms_rhs);
	i_trsm( side, trans_tri_rhs, alpha, tri_rhs, gms_lhs );
}

}	// end namespace

// gm_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) (left) (BLAS xTRMM).
void LinAlgPack::M_StMtM(GenMatrix* gm_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(  tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1
					 , gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2				);
	gm_lhs->resize(	  rows(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1)
					, cols(gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2)				);
	i_trmm_alt(BLAS_Cpp::left,alpha,tri_rhs1,trans_rhs1,gms_rhs2,trans_rhs2,&(*gm_lhs)());
}

// gms_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) (left) (BLAS xTRMM).
void LinAlgPack::M_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1
							, gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2				);
	i_trmm_alt(BLAS_Cpp::left,alpha,tri_rhs1,trans_rhs1,gms_rhs2,trans_rhs2,gms_lhs);
}

// gm_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM).
void LinAlgPack::M_StMtM(GenMatrix* gm_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(  gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
					 , tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2		);
	gm_lhs->resize(	  rows(gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1)
					, cols(tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2)		);
	i_trmm_alt(BLAS_Cpp::right,alpha,tri_rhs2,trans_rhs2,gms_rhs1,trans_rhs1,&(*gm_lhs)());
}

// gms_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM).
void LinAlgPack::M_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
							, tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2 );
	i_trmm_alt(BLAS_Cpp::right,alpha,tri_rhs2,trans_rhs2,gms_rhs1,trans_rhs1,gms_lhs);
}

// gms_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) + beta * gms_lhs (left).
void LinAlgPack::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	// If op(gms_rhs2) == gms_lhs and beta = 0.0 then this is a direct call to the BLAS.
	if(	   gms_lhs->overlap(gms_rhs2) == SAME_MEM && trans_rhs2 == BLAS_Cpp::no_trans
		&& beta == 0.0 )
	{
		i_trmm(BLAS_Cpp::left,trans_rhs1,alpha,tri_rhs1,gms_lhs);
	}
	else {
		// This is something else so the alias problem is not handled.
		if(beta != 1.0) Mt_S(gms_lhs,beta);
		GenMatrix tmp;
		M_StMtM(&tmp,alpha,tri_rhs1,trans_rhs1,gms_rhs2,trans_rhs2);
		Mp_StM(gms_lhs,1.0,tmp,BLAS_Cpp::no_trans);
	}
}

// gms_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) + beta * gms_lhs (right).
void LinAlgPack::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	// If op(gms_rhs1) == gms_lhs and beta = 0.0 then this is a direct call to the BLAS.
	if(	   gms_lhs->overlap(gms_rhs1) == SAME_MEM && trans_rhs1 == BLAS_Cpp::no_trans
		&& beta == 0.0 )
	{
		i_trmm(BLAS_Cpp::right,trans_rhs2,alpha,tri_rhs2,gms_lhs);
	}
	else {
		// This is something else so the alias problem is not handled.
		if(beta != 1.0) Mt_S(gms_lhs,beta);
		GenMatrix tmp;
		M_StMtM(&tmp,alpha,gms_rhs1,trans_rhs1,tri_rhs2,trans_rhs2);
		Mp_StM(gms_lhs,1.0,tmp,BLAS_Cpp::no_trans);
	}
}

// gm_lhs = alpha * inv(op(tri_rhs1)) * op(gms_rhs2) (left) (BLAS xTRSM).
void LinAlgPack::M_StInvMtM(GenMatrix* gm_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(  tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1
					 , gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2					);
	gm_lhs->resize(	  rows(tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1)
					, cols(gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2)				);
	i_trsm_alt(BLAS_Cpp::left,alpha,tri_rhs1,trans_rhs1,gms_rhs2,trans_rhs2,&(*gm_lhs)());
}

// gms_lhs = alpha * inv(op(tri_rhs1)) * op(gms_rhs2) (left) (BLAS xTRSM).
void LinAlgPack::M_StInvMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, tri_rhs1.gms().rows(), tri_rhs1.gms().cols(), trans_rhs1
							, gms_rhs2.rows(), gms_rhs2.cols(), trans_rhs2 );
	i_trsm_alt(BLAS_Cpp::left,alpha,tri_rhs1,trans_rhs1,gms_rhs2,trans_rhs2,gms_lhs);
}

// gm_lhs = alpha * op(gms_rhs1) * inv(op(tri_rhs2)) (right) (BLAS xTRSM).
void LinAlgPack::M_StMtInvM(GenMatrix* gm_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
						, tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2		);
	gm_lhs->resize(	  rows(gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1)
					, cols(tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2)		);
	i_trsm_alt(BLAS_Cpp::right,alpha,tri_rhs2,trans_rhs2,gms_rhs1,trans_rhs1,&(*gm_lhs)());
}

// gms_lhs = alpha * op(gms_rhs1) * inv(op(tri_rhs2)) (right) (BLAS xTRSM).
void LinAlgPack::M_StMtInvM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, gms_rhs1.rows(), gms_rhs1.cols(), trans_rhs1
							, tri_rhs2.gms().rows(), tri_rhs2.gms().cols(), trans_rhs2 );
	i_trsm_alt(BLAS_Cpp::right,alpha,tri_rhs2,trans_rhs2,gms_rhs1,trans_rhs1,gms_lhs);
}
