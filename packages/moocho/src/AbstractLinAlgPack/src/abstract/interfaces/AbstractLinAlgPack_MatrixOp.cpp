// //////////////////////////////////////////////////////////////
// MatrixWithOp.cpp
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

#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpSubView.h"
#include "AbstractLinAlgPack/include/MatrixPermAggr.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/Permutation.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/SpVectorView.h"
#include "AbstractLinAlgPack/include/EtaVector.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

void MatrixWithOp::zero_out()
{
	THROW_EXCEPTION(
		true, std::logic_error, "MatrixWithOp::zero_out(): "
		"Error, this method as not been defined by the subclass \'"
		<<typeid(*this).name()<<"\'" );
}

void MatrixWithOp::Mt_S(value_type alpha)
{
	THROW_EXCEPTION(
		true, std::logic_error, "MatrixWithOp::Mt_S(): "
		"Error, this method as not been defined by the subclass \'"
		<<typeid(*this).name()<<"\'" );
}

MatrixWithOp& MatrixWithOp::operator=(const MatrixWithOp& M)
{
	const bool assign_to_self = dynamic_cast<const void*>(this) == dynamic_cast<const void*>(&M);
	THROW_EXCEPTION(
		!assign_to_self, std::logic_error
		,"MatrixWithOp::operator=(M) : Error, this is not assignment "
		"to self and this method is not overridden for the subclass \'"
		<< typeid(*this).name() << "\'" );
	return *this; // assignment to self
}

std::ostream& MatrixWithOp::output(std::ostream& out) const
{
	const size_type m = this->rows(), n = this->cols();
	VectorSpace::vec_mut_ptr_t
		row_vec = space_rows().create_member(); // dim() == n
	out << m << " " << n << std::endl;
	for( size_type i = 1; i <= m; ++i ) {
		LinAlgOpPack::V_StMtV( row_vec.get(), 1.0, *this, BLAS_Cpp::trans, EtaVector(i,m)() );
		row_vec->output(out,false,true);
	}
	return out;
}

// Clone

MatrixWithOp::mat_mut_ptr_t
MatrixWithOp::clone()
{
	return MemMngPack::null;
}

MatrixWithOp::mat_ptr_t
MatrixWithOp::clone() const
{
	return MemMngPack::null;
}

// Subview

MatrixWithOp::mat_ptr_t
MatrixWithOp::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	namespace rcp = MemMngPack;

	if( 
		( ( row_rng.lbound() == 1 && row_rng.ubound() == this->rows() )
		  || row_rng.full_range() )
		&&
		( ( col_rng.lbound() == 1 && col_rng.ubound() == this->cols() )
		  || row_rng.full_range() )
		) 
	{
		return rcp::rcp(this,false); // don't clean up memory
	}
	return rcp::rcp(
		new MatrixWithOpSubView(
			rcp::rcp(const_cast<MatrixWithOp*>(this),false) // don't clean up memory
			,row_rng,col_rng ) );
}

// Permuted views

MatrixWithOp::mat_ptr_t
MatrixWithOp::perm_view(
	const Permutation          *P_row
	,const index_type          row_part[]
	,int                       num_row_part
	,const Permutation         *P_col
	,const index_type          col_part[]
	,int                       num_col_part
	) const
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new MatrixPermAggr(
			rcp::rcp(this,false)
			,rcp::rcp(P_row,false)
			,rcp::rcp(P_col,false)
			,rcp::null
			) );
}

MatrixWithOp::mat_ptr_t
MatrixWithOp::perm_view_update(
	const Permutation          *P_row
	,const index_type          row_part[]
	,int                       num_row_part
	,const Permutation         *P_col
	,const index_type          col_part[]
	,int                       num_col_part
	,const mat_ptr_t           &perm_view
	) const
{
	return this->perm_view(
		P_row,row_part,num_row_part
		,P_col,col_part,num_col_part );
}

// Level-1 BLAS

bool MatrixWithOp::Mp_StM(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;

	// Give m_lhs a chance to implement this method
	if(m_lhs->Mp_StM(alpha,*this,trans_rhs))
		return true;

	// We msut try to implement the method
	MultiVectorMutable
		*m_mut_lhs = dynamic_cast<MultiVectorMutable*>(m_lhs);
	if(!m_mut_lhs)
		return false; // Could not implement the method :-(
		
#ifdef _DEBUG
	THROW_EXCEPTION(
		!m_lhs->space_rows().is_compatible(
			trans_rhs == no_trans ? this->space_rows() : this->space_cols() )
		|| !m_lhs->space_cols().is_compatible(
			trans_rhs == no_trans ? this->space_cols() : this->space_rows() )
		, IncompatibleMatrices
		,"MatrixWithOp::Mp_StM(m_lhs,...): Error, m_lhs of type \'"<<typeid(*m_lhs).name()<<"\' "
		<<"is not compatible with this of type \'"<<typeid(*m_lhs).name()<<"\'" );
#endif

	const size_type
		rows = BLAS_Cpp::rows( m_lhs->rows(), m_lhs->cols(), trans_rhs ),
		cols = BLAS_Cpp::cols( m_lhs->rows(), m_lhs->cols(), trans_rhs );
	for( size_type j = 1; j <= cols; ++j )
		AbstractLinAlgPack::Vp_StMtV( m_mut_lhs->col(j).get(), alpha, *this, trans_rhs, EtaVector(j,cols)() );
	// ToDo: consider row access?
	return true;
}

bool MatrixWithOp::Mp_StM(
	value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
	return false;
}

bool MatrixWithOp::Mp_StMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	// Give m_lhs a chance to implement this method
	if(m_lhs->Mp_StMtP(alpha,*this,M_trans,P_rhs,P_rhs_trans))
		return true;
	assert(0); // ToDo: Implement default!
	return false;
}

bool MatrixWithOp::Mp_StMtP(
	value_type alpha
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	return false;
}

bool MatrixWithOp::Mp_StPtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	// Give m_lhs a chance to implement this method
	if(m_lhs->Mp_StPtM(alpha,P_rhs,P_rhs_trans,*this,M_trans))
		return true;
	assert(0); // ToDo: Implement default!
	return false;
}

bool MatrixWithOp::Mp_StPtM(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	return false;
}

bool MatrixWithOp::Mp_StPtMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	// Give m_lhs a chance to implement this method
	if(m_lhs->Mp_StPtMtP(alpha,P_rhs1,P_rhs1_trans,*this,M_trans,P_rhs2,P_rhs2_trans))
		return true;
	assert(0); // ToDo: Implement default!
	return false;
}

bool MatrixWithOp::Mp_StPtMtP(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	return false;
}

// Level-2 BLAS

void MatrixWithOp::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	Vp_MtV_assert_compatibility(v_lhs,*this,trans_rhs1,sv_rhs2 );
	if( sv_rhs2.nz() ) {
		VectorSpace::vec_mut_ptr_t
			v_rhs2 = (trans_rhs1 == BLAS_Cpp::no_trans
					  ? this->space_rows()
					  : this->space_cols()
				).create_member();
		v_rhs2->set_sub_vector(sub_vec_view(sv_rhs2));
		this->Vp_StMtV(v_lhs,alpha,trans_rhs1,*v_rhs2,beta);
	}
	else {
		Vt_S( v_lhs, beta );
	}
}

void MatrixWithOp::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& vs_rhs3, value_type beta) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	assert(0); // ToDo: Implement!
}

value_type MatrixWithOp::transVtMtV(
	const VectorWithOp& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorWithOp& vs_rhs3) const
{
	assert(0); // ToDo: Implement!
	return 0.0;
}

value_type MatrixWithOp::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	assert(0); // ToDo: Implement!
	return 0.0;
}

void MatrixWithOp::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert(0); // ToDo: Implement!
}

// Level-3 BLAS

bool MatrixWithOp::Mp_StMtM(
	MatrixWithOp* C, value_type a
	, BLAS_Cpp::Transp A_trans, const MatrixWithOp& B
	, BLAS_Cpp::Transp B_trans, value_type b) const
{
	if(B.Mp_StMtM(C,a,*this,A_trans,B_trans,b))
		return true;
	else if(C->Mp_StMtM(a,*this,A_trans,B,B_trans,b))
		return true;
	//
	// C = b*C + a*op(A)*op(B), where A == *this
	//
	// We will perform this by column as:
	//
	//   C(:,j) = b*C(:,j) + a*op(A)*op(B)*e(j), for j = 1...C.cols()
	//
	//   by performing:
	//
	//        t = op(B)*e(j)
	//        C(:,j) = b*C(:,j) + a*op(A)*t
	//
	Mp_MtM_assert_compatibility(C,BLAS_Cpp::no_trans,*this,A_trans,B,B_trans);
	MultiVectorMutable *Cmv = dynamic_cast<MultiVectorMutable*>(C);
	if(!Cmv)
		return false; // MultiVectorMutable is not supported!
	if( !(Cmv->access_by() & MultiVector::COL_ACCESS) )
		return false; // Column access is not supported!
	VectorSpace::vec_mut_ptr_t
		t = ( B_trans == BLAS_Cpp::no_trans ? B.space_cols() : B.space_rows() ).create_member();
	const index_type
		C_rows = Cmv->rows(),
		C_cols = Cmv->cols();
	for( index_type j = 1; j <= C_cols; ++j ) {
		// t = op(B)*e(j)
		LinAlgOpPack::V_MtV( t.get(), B, B_trans, EtaVector(j,C_cols)	);	
		// C(:,j) = a*op(A)*t + b*C(:,j)
		this->Vp_StMtV( Cmv->col(j).get(), a, A_trans, *t, b );
	}
	return false;
}

bool MatrixWithOp::Mp_StMtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	return false;
}

bool MatrixWithOp::Mp_StMtM(
	value_type alpha
	,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
	,value_type beta )
{
	return false;
}

void MatrixWithOp::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert(0); // ToDo: Implement!
}

// overridden from MatrixBase

size_type MatrixWithOp::rows() const
{
	return this->space_cols().dim();
}

size_type MatrixWithOp::cols() const
{
	return this->space_rows().dim();
}

} // end namespace AbstractLinAlgPack

// Non-member functions

// level-1 BLAS

void AbstractLinAlgPack::Mp_StM(
	MatrixWithOp* mwo_lhs, value_type alpha, const MatrixWithOp& M_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	const bool success = M_rhs.Mp_StM(mwo_lhs,alpha,trans_rhs);
	THROW_EXCEPTION(
		!success, MatrixWithOp::MethodNotImplemented
		,"MatrixWithOp::Mp_StM(...) : Error, mwo_lhs of type \'"
		<< typeid(*mwo_lhs).name() << "\' does not support the "
		"\'MultiVectorMutable\' interface!" );
}

void AbstractLinAlgPack::Mp_StMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	const bool success = M_rhs.Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans);
	THROW_EXCEPTION(
		!success, MatrixWithOp::MethodNotImplemented
		,"MatrixWithOp::Mp_StMtP(...) : Error, mwo_lhs of type \'"
		<< typeid(*mwo_lhs).name() << "\' does not support the "
		"\'MultiVectorMutable\' interface!" );
}

void AbstractLinAlgPack::Mp_StPtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	const bool success = M_rhs.Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans);
	THROW_EXCEPTION(
		!success, MatrixWithOp::MethodNotImplemented
		,"MatrixWithOp::Mp_StPtM(...) : Error, mwo_lhs of type \'"
		<< typeid(*mwo_lhs).name() << "\' does not support the "
		"\'MultiVectorMutable\' interface!" );
}

void AbstractLinAlgPack::Mp_StPtMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	const bool success = M_rhs.Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,trans_rhs,P_rhs2,P_rhs2_trans);
	THROW_EXCEPTION(
		!success, MatrixWithOp::MethodNotImplemented
		,"MatrixWithOp::Mp_StPtMtP(...) : Error, mwo_lhs of type \'"
		<< typeid(*mwo_lhs).name() << "\' does not support the "
		"\'MultiVectorMutable\' interface!" );
}

// level-3 blas

void AbstractLinAlgPack::Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta )
{
	const bool success = mwo_rhs1.Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
	THROW_EXCEPTION(
		!success, MatrixWithOp::MethodNotImplemented
		,"MatrixWithOp::Mp_StMtM(...) : Error, mwo_lhs of type \'"
		<< typeid(*mwo_lhs).name() << "\' does not support the "
		"\'MultiVectorMutable\' interface or does not support column access!" );
}
