// ////////////////////////////////////////////////////////////////////////
// MatrixWithOpSerial.cpp
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

// ToDo: 7/27/99: Give these default implementations and test them.

#include <typeinfo>

#include "SparseLinAlgPack/include/MatrixWithOpSerial.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "SparseLinAlgPack/include/VectorWithOpGetSparse.h"
#include "SparseLinAlgPack/include/MatrixWithOpGetGMSMutable.h"
#include "SparseLinAlgPack/include/MatrixWithOpGetGMSTri.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpGetGMSSymMutable.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "AbstractLinAlgPack/include/EtaVector.h"
#include "AbstractLinAlgPack/include/GenPermMatrixSlice.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "WorkspacePack.h"
#include "dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Mp_StM;
}

namespace SparseLinAlgPack {

// Level-1 BLAS

void MatrixWithOpSerial::Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	LinAlgPack::Mp_M_assert_sizes( gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), trans_rhs );
	const size_type
		m = gms_lhs->rows(),
		n = gms_lhs->cols();
	//
	// Use sparse matrix-vector multiplication to perform this operation.
	// C += a * B = a * B * I = [ a*B*e(1), a*B*e(2), ..., a*B*e(m) ]
	//
	SpVector rhs;
	rhs.uninitialized_resize( n, 1, 1 );
	for( size_type j = 1; j <=n; ++j ) {
		rhs.begin()->initialize( j, 1.0 );	// e(j)
		this->Vp_StMtV( &gms_lhs->col(j), alpha, trans_rhs, rhs(), 1.0 );
	}
}

void MatrixWithOpSerial::Mp_StMtP(GenMatrixSlice* C, value_type a
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	) const 
{
	// C += a * op(M) * op(P)
	assert(0);	// Implement this!
}

void MatrixWithOpSerial::Mp_StPtM(GenMatrixSlice* C, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	) const 
{
	// C += a * op(P) * op(M)
	assert(0);	// Implement this!
}

void MatrixWithOpSerial::Mp_StPtMtP( GenMatrixSlice* C, value_type a
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	) const
{
	// C += a * op(P1) * op(M) * op(P2)
	assert(0);	// Implement this!
}

// Level-2 BLAS

void MatrixWithOpSerial::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	Vp_MtV_assert_sizes( vs_lhs->dim(), rows(), cols(), trans_rhs1, sv_rhs2.dim() );
	if( !sv_rhs2.nz() ) {
		// vs_lhs = beta * vs_lhs
		if(beta==0.0)      *vs_lhs = 0.0;
		else if(beta!=1.0) LinAlgPack::Vt_S(vs_lhs,beta);
	}
	else {
		// Convert to dense by default.
		if( sv_rhs2.dim() == sv_rhs2.nz() && sv_rhs2.is_sorted() ) {
			const VectorSlice vs_rhs2 = SparseLinAlgPack::dense_view(sv_rhs2);
			this->Vp_StMtV( vs_lhs, alpha, trans_rhs1, vs_rhs2, beta );
		}
		else {
			Vector v_rhs2;
			LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
			this->Vp_StMtV( vs_lhs, alpha, trans_rhs1, v_rhs2(), beta );
		}
	}
}

void MatrixWithOpSerial::Vp_StPtMtV(VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const VectorSlice& x, value_type b) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	wsp::Workspace<value_type> t_ws(wss,BLAS_Cpp::cols(P.rows(),P.cols(),P_trans));
	VectorSlice                t(&t_ws[0],t_ws.size());
    LinAlgOpPack::V_StMtV(&t,a,*this,M_trans,x);
	LinAlgOpPack::Vp_MtV( y, P, P_trans, t, b ); 
}

void MatrixWithOpSerial::Vp_StPtMtV(VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	wsp::Workspace<value_type> t_ws(wss,BLAS_Cpp::cols(P.rows(),P.cols(),P_trans));
	VectorSlice                t(&t_ws[0],t_ws.size());
    LinAlgOpPack::V_StMtV(&t,a,*this,M_trans,x);
	LinAlgOpPack::Vp_MtV( y, P, P_trans, t, b ); 
}

value_type MatrixWithOpSerial::transVtMtV(const VectorSlice& x1
	, BLAS_Cpp::Transp M_trans, const VectorSlice& x2) const
{
	LinAlgPack::Vp_MtV_assert_sizes( x1.dim(), rows(), cols(), M_trans, x2.dim() );
	Vector tmp(x1.dim());
	this->Vp_StMtV( &tmp(), 1.0, M_trans, x2, 0.0 );
	return LinAlgPack::dot( x1, tmp() );
}

value_type MatrixWithOpSerial::transVtMtV(const SpVectorSlice& x1
	, BLAS_Cpp::Transp M_trans, const SpVectorSlice& x2) const
{
	LinAlgPack::Vp_MtV_assert_sizes( x1.dim(), rows(), cols(), M_trans, x2.dim() );
	if( !x1.nz() || !x2.nz() ) {
		return 0.0;
	}
	else {
		if( x1.overlap(x2) == LinAlgPack::SAME_MEM && x1.dim() == x1.nz() && x1.is_sorted()  ) {
			const VectorSlice x1_d = SparseLinAlgPack::dense_view(x1);
			return this->transVtMtV( x1_d, M_trans, x1_d );
		}
		Vector tmp(x1.dim());
		this->Vp_StMtV( &tmp(), 1.0, M_trans, x2, 0.0 );
		return SparseLinAlgPack::dot( x1, tmp() );
	}
}

void MatrixWithOpSerial::syr2k(
	 BLAS_Cpp::Transp M_trans_in, value_type a
	, const GenPermMatrixSlice& P1_in, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2_in, BLAS_Cpp::Transp P2_trans
	, value_type b, sym_gms* S ) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	//
	// S = b * S
	//
	// S += a*op(P1')*op(M)*op(P2) + a*op(P2')*op(M')*op(P1)
	//
	// We will start by renaming P1 and P2 such that op(P1).rows() >= op(P2).rows().
	// This is because we are going to store some temparary vectors and we don't
	// want them to be too big.
	//
	// We will perform the above operation by working with columns of:
	//
	//    op(P1)(:,j(k)) = e(i(k)) <: R^n
	//
	// Then for each column in op(P1) we will perform:
	//
	//
	// for k = 1...P1.nz()
	//
	//              [    .     ]
	//     S += a * [ e(i(k))' ] * op(M)*op(P2) + a * op(P2') * op(M') * [  ...  e(i(k))  ...  ] 
	//              [    .     ]
	//                row j(k)                                                    col j(k)
	//     =>
	//              [  .   ]
	//     S += a * [ y_k' ] + a * [  ...  y_k  ...  ] 
	//              [  .   ]
	//               row j(k)            col j(k)
	//
	//     where: y_k = a * op(P2') * op(M') * e(i(k)) <: R^m
	// 
	// Of course above we only need to set the row and column elements for S(j(k),:) and S(:,j(k))
	// for the portion of the symmetric S that is being stored.
	//
	const size_type
		M_rows  = this->rows(),
		M_cols  = this->cols(),
		P1_cols = cols( P1_in.rows(), P1_in.cols(), P1_trans );
	LinAlgPack::MtM_assert_sizes(
		M_rows, M_cols, trans_not(M_trans_in)
		, P1_in.rows(), P1_in.cols(), P1_trans );
	LinAlgPack::MtM_assert_sizes(
		M_rows, M_cols, M_trans_in
		, P2_in.rows(), P2_in.cols(), P2_trans );
	LinAlgPack::Mp_M_assert_sizes(
		S->rows(), S->cols(), no_trans
		, P1_cols, P1_cols, no_trans );
	// Rename P1 and P2 so that op(P1).rows() >= op(P2).rows()
	const bool
		reorder_P1_P2 = ( rows( P1_in.rows(), P1_in.cols(), P1_trans ) 
						  < rows( P2_in.rows(), P2_in.cols(), P2_trans ) );
	const GenPermMatrixSlice
		&P1 = reorder_P1_P2 ? P2_in : P1_in,
		&P2 = reorder_P1_P2 ? P1_in : P2_in;
	const BLAS_Cpp::Transp
		M_trans  = reorder_P1_P2 ? trans_not(M_trans_in) : M_trans_in;
	// Set rows and columns of S
	const size_type
		m  = S->rows(),
		n = rows( P1.rows(), P1.cols(), P1_trans );
	Vector y_k_store(m); // ToDo: use workspace allocator!
	VectorSlice y_k = y_k_store();
	for( GenPermMatrixSlice::const_iterator P1_itr = P1.begin(); P1_itr != P1.end(); ++P1_itr )
	{
		const size_type
			i_k = P1_trans == no_trans ? P1_itr->row_i() : P1_itr->col_j(),
			j_k = P1_trans == no_trans ? P1_itr->col_j() : P1_itr->row_i();
		// e(i(k))
		EtaVector
			e_i_k(i_k,n);
		// y_k = op(P2')*op(M')*e(i(k))
		SparseLinAlgPack::Vp_StPtMtV( &y_k, 1.0, P2, trans_not(P2_trans), *this, trans_not(M_trans), e_i_k(), 0.0 );
		// S(j(k),:) += a*y_k'
		if( S->uplo() == BLAS_Cpp::upper )
			LinAlgPack::Vp_StV( &S->gms().row(j_k)(1,j_k), a, y_k(1,j_k) );
		else
			LinAlgPack::Vp_StV( &S->gms().row(j_k)(j_k,m), a, y_k(j_k,m) );
		// S(:,j(k)) += a*y_k
		if( S->uplo() == BLAS_Cpp::upper )
			LinAlgPack::Vp_StV( &S->gms().col(j_k)(1,j_k), a, y_k(1,j_k) );
		else
			LinAlgPack::Vp_StV( &S->gms().col(j_k)(j_k,m), a, y_k(j_k,m) );
	}
}

// Level-3 BLAS

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* C, value_type a
	, BLAS_Cpp::Transp A_trans, const GenMatrixSlice& B
	, BLAS_Cpp::Transp B_trans, value_type b) const
{
	LinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	//
	// C = b*C + a*op(A)*op(B)
	// 
	// C.col(j) = b*col(j) + a*op(A)*op(B).col(j)
	//

	// ToDo: Add the ability to also perform by row if faster

	for( size_type j = 1; j <= C->cols(); ++j )
		SparseLinAlgPack::Vp_StMtV( &C->col(j), a, *this, A_trans, LinAlgPack::col( B, B_trans, j ), b );
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* C, value_type a, const GenMatrixSlice& A
	, BLAS_Cpp::Transp A_trans, BLAS_Cpp::Transp B_trans, value_type b) const
{
	LinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
		, A.rows(), A.cols(), A_trans, rows(), cols(), B_trans );
	//
	// C = b*C + a*op(A)*op(B)
	// 
	// C.row(i) = b*row(i) + a*op(B)'*op(A).row(i)
	//

	// ToDo: Add the ability to also perform by column if faster

	for( size_type i = 1; i <= C->rows(); ++i )
		SparseLinAlgPack::Vp_StMtV( &C->row(i), a, *this, BLAS_Cpp::trans_not(A_trans)
			, LinAlgPack::row(A,A_trans,i) , b );
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* C, value_type a
	, BLAS_Cpp::Transp A_trans, const MatrixWithOpSerial& B
	, BLAS_Cpp::Transp B_trans, value_type b) const
{
	using LinAlgOpPack::assign;
	// C = b*C + a*op(A)*op(B)
	LinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
	// Convert one of the matrices to dense, which ever one is the smallest!
	const size_type
		size_A = rows() * cols(),
		size_B = B.rows() * B.cols();
	if( size_A < size_B ) {
		GenMatrix A_dense;
		assign( &A_dense, *this, BLAS_Cpp::no_trans );
		SparseLinAlgPack::Mp_StMtM( C, a, A_dense(), A_trans, B, B_trans, b );
	}
	else {
		GenMatrix B_dense;
		assign( &B_dense, B, BLAS_Cpp::no_trans );
		SparseLinAlgPack::Mp_StMtM( C, a, *this, A_trans, B_dense(), B_trans, b );
	}
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const sym_gms& sym_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert(0); // Todo: Implement!
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const sym_gms& sym_rhs1
	, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert(0); // Todo: Implement!
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert(0); // Todo: Implement!
}

void MatrixWithOpSerial::Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert(0); // Todo: Implement!
}

void MatrixWithOpSerial:: syrk(
	BLAS_Cpp::Transp M_trans, value_type a
	, value_type b, sym_gms* S ) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	//
	// S = b*S + a*op(M)*op(M')
	//
	const size_type
		M_rows = this->rows(),
		M_cols = this->cols(),
		opM_rows = rows( M_rows, M_cols, M_trans ),
		opM_cols = cols( M_rows, M_cols, M_trans ),
		m = opM_rows;
	LinAlgPack::Mp_MtM_assert_sizes(
		S->rows(), S->cols(), no_trans
		,M_rows, M_cols, M_trans
		,M_rows, M_cols, trans_not(M_trans) );
	// S = b*S
	LinAlgPack::Mt_S( &LinAlgPack::nonconst_tri_ele(S->gms(),S->uplo()), b );
	//
	// Form this matrix one column at a time by multiplying by e(j):
	//
	// S(:,j) += a*op(M)*(op(M')*e(j))
	//
	//    j = 1 ... opM_rows
	//
	wsp::Workspace<value_type> t1_ws(wss,opM_cols),
		                       t2_ws(wss,opM_rows);
	VectorSlice                t1(&t1_ws[0],t1_ws.size()),
		                       t2(&t2_ws[0],t2_ws.size());
	for( size_type j = 1; j <= opM_rows; ++j ) {
		EtaVector e_j(j,opM_rows);
		LinAlgOpPack::V_MtV(&t1,*this,trans_not(M_trans),e_j()); // t1 = op(M')*e(j)
		LinAlgOpPack::V_MtV(&t2,*this,M_trans,t1);               // t2 = op(M)*t1
		// S(j,:) += a*t2' 
		if( S->uplo() == BLAS_Cpp::upper )
			LinAlgPack::Vp_StV( &S->gms().row(j)(1,j), a, t2(1,j) );
		else
			LinAlgPack::Vp_StV( &S->gms().row(j)(j,m), a, t2(j,m) );
		// S(:,j) += a*t2
		if( S->uplo() == BLAS_Cpp::upper )
			LinAlgPack::Vp_StV( &S->gms().col(j)(1,j), a, t2(1,j) );
		else
			LinAlgPack::Vp_StV( &S->gms().col(j)(j,m), a, t2(j,m) );
	}
}

// Overridden from MatrixWithOp
	
const VectorSpace& MatrixWithOpSerial::space_cols() const
{
	const size_type rows = this->rows();
	if(space_cols_.dim() != rows)
		space_cols_.initialize(rows);
	return space_cols_;
}
	
const VectorSpace& MatrixWithOpSerial::space_rows() const
{
	const size_type cols = this->cols();
	if(space_rows_.dim() != cols)
		space_rows_.initialize(cols);
	return space_rows_;
}
	
std::ostream& MatrixWithOpSerial::output(std::ostream& out) const {
	GenMatrix tmp( 0.0, rows(), cols() );
	this->Mp_StM( &tmp(), 1.0 , BLAS_Cpp::no_trans );
	return out << tmp;
}

bool MatrixWithOpSerial::Mp_StM(
	MatrixWithOp* mwo_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs
	) const
{
	MatrixWithOpGetGMSMutable
		*mwo_gms_lhs = dynamic_cast<MatrixWithOpGetGMSMutable*>(mwo_lhs);
	if(!mwo_gms_lhs)
		return MatrixWithOp::Mp_StM(mwo_lhs,alpha,trans_rhs); // boot it!
	this->Mp_StM( &MatrixDenseMutableEncap(mwo_gms_lhs)(), alpha, trans_rhs );
	return true;
}
	
bool MatrixWithOpSerial::Mp_StMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	MatrixWithOpGetGMSMutable
		*mwo_gms_lhs = dynamic_cast<MatrixWithOpGetGMSMutable*>(mwo_lhs);
	if(!mwo_gms_lhs)
		return MatrixWithOp::Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans); // boot it!
	this->Mp_StMtP(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,M_trans,P_rhs,P_rhs_trans);
	return true;
}
	
bool MatrixWithOpSerial::Mp_StPtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	MatrixWithOpGetGMSMutable
		*mwo_gms_lhs = dynamic_cast<MatrixWithOpGetGMSMutable*>(mwo_lhs);
	if(!mwo_gms_lhs)
		return MatrixWithOp::Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans); // boot it!
	this->Mp_StPtM(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,P_rhs,P_rhs_trans,M_trans);
	return true;
}
	
bool MatrixWithOpSerial::Mp_StPtMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	MatrixWithOpGetGMSMutable
		*mwo_gms_lhs = dynamic_cast<MatrixWithOpGetGMSMutable*>(mwo_lhs);
	if(!mwo_gms_lhs)
		return MatrixWithOp::Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans); // boot it!
	this->Mp_StPtMtP(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans);
	return true;
}
	
void MatrixWithOpSerial::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorWithOp& v_rhs2, value_type beta) const
{
	VectorDenseMutableEncap       vs_lhs(*v_lhs);
	const VectorWithOpGetSparse   *sv_rhs2 = dynamic_cast<const VectorWithOpGetSparse*>(&v_rhs2);
	if(sv_rhs2)
		this->Vp_StMtV( &vs_lhs(), alpha, trans_rhs1, VectorSparseEncap(*sv_rhs2)(), beta );
	VectorDenseEncap              vs_rhs2(v_rhs2);
	this->Vp_StMtV( &vs_lhs(), alpha, trans_rhs1, vs_rhs2(), beta );	
}

void MatrixWithOpSerial::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	VectorDenseMutableEncap       vs_lhs(*v_lhs);
	this->Vp_StMtV( &vs_lhs(), alpha, trans_rhs1, sv_rhs2, beta );	
}
	
void MatrixWithOpSerial::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& v_rhs3, value_type beta) const
{
	VectorDenseMutableEncap       vs_lhs(*v_lhs);
	const VectorWithOpGetSparse   *sv_rhs3 = dynamic_cast<const VectorWithOpGetSparse*>(&v_rhs3);
	if(sv_rhs3) {
		this->Vp_StPtMtV( &vs_lhs(), alpha, P_rhs1, P_rhs1_trans, M_rhs2_trans, VectorSparseEncap(*sv_rhs3)(), beta );
		return;
	}
	VectorDenseEncap              vs_rhs3(v_rhs3);
	this->Vp_StPtMtV( &vs_lhs(), alpha, P_rhs1, P_rhs1_trans, M_rhs2_trans, vs_rhs3(), beta );	
}
	
void MatrixWithOpSerial::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	VectorDenseMutableEncap       vs_lhs(*v_lhs);
	this->Vp_StPtMtV( &vs_lhs(), alpha, P_rhs1, P_rhs1_trans, M_rhs2_trans, sv_rhs3, beta );	
}
	
value_type MatrixWithOpSerial::transVtMtV(
	const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorWithOp& v_rhs3) const
{
	VectorDenseEncap              vs_rhs1(v_rhs1);
	VectorDenseEncap              vs_rhs3(v_rhs3);
	return this->transVtMtV(vs_rhs1(),trans_rhs2,vs_rhs3());
}
	
void MatrixWithOpSerial::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymWithOp* symwo_lhs ) const
{
	MatrixSymWithOpGetGMSSymMutable
		*symwo_gms_lhs = dynamic_cast<MatrixSymWithOpGetGMSSymMutable*>(symwo_lhs);
	if(!symwo_gms_lhs) {
		MatrixWithOp::syr2k(M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs); // Boot it
		return;
	}
	this->syr2k(
		M_trans,alpha,P1,P1_trans,P2,P2_trans,beta
		,&MatrixDenseSymMutableEncap(symwo_gms_lhs)()
		 );
}

bool MatrixWithOpSerial::Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	MatrixWithOpGetGMSMutable
		*mwo_gms_lhs = dynamic_cast<MatrixWithOpGetGMSMutable*>(mwo_lhs);
	if(mwo_gms_lhs) {
		// Try to match the rhs arguments to known serial interfaces
		if(const MatrixWithOpGetGMS* mwo_gms_rhs2 = dynamic_cast<const MatrixWithOpGetGMS*>(&mwo_rhs2)) {
			this->Mp_StMtM(
				&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
				,MatrixDenseEncap(*mwo_gms_rhs2)(),trans_rhs2,beta );
			return true;
		}
		if(const MatrixSymWithOpGetGMSSym* mwo_sym_gms_rhs2 = dynamic_cast<const MatrixSymWithOpGetGMSSym*>(&mwo_rhs2)) {
			this->Mp_StMtM(
				&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
				,MatrixDenseEncap(*mwo_sym_gms_rhs2)(),trans_rhs2,beta );
			return true;
		}
		if(const MatrixWithOpGetGMSTri* mwo_tri_gms_rhs2 = dynamic_cast<const MatrixWithOpGetGMSTri*>(&mwo_rhs2)) {
			this->Mp_StMtM(
				&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
				,MatrixDenseEncap(*mwo_tri_gms_rhs2)(),trans_rhs2,beta );
			return true;
		}
		// If we get here, the matrix arguments did not match up so we have to give up (I think?)
	}
	// Let the default implementation try to find matrix arguments that can handle this!
	return MatrixWithOp::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // Boot it!
}

bool MatrixWithOpSerial::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* symwo_lhs ) const
{
	MatrixSymWithOpGetGMSSymMutable
		*symwo_gms_lhs = dynamic_cast<MatrixSymWithOpGetGMSSymMutable*>(symwo_lhs);
	if(!symwo_gms_lhs) {
		return MatrixWithOp::syrk(M_trans,alpha,beta,symwo_lhs); // Boot it
	}
	this->syrk(M_trans,alpha,beta,&MatrixDenseSymMutableEncap(symwo_gms_lhs)());
	return true;
}

}	// end namespace SparseLinAlgPack
