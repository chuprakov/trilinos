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

#include <iostream> // Debuggin only

#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

MatrixSymDiagonalStd::MatrixSymDiagonalStd(
	const VectorSpace::vec_mut_ptr_t& diag
	,bool                             unique
	)
{
	this->initialize(diag,unique);
//	std::cerr << "MatrixSymDiagonalStd::rows() = " << this->rows() << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagonalStd::nz() = "   << this->nz()   << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagonalStd::cols() = " << this->cols() << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagonalStd::nz() = "   << this->nz()   << std::endl; // Debugging
}

void MatrixSymDiagonalStd::initialize(
	const VectorSpace::vec_mut_ptr_t& diag
	,bool                             unique
	)
{
	diag_   = diag;   // lazy copy!
	unique_ = unique;
}

VectorWithOpMutable& MatrixSymDiagonalStd::diag()
{
	copy_unique();
	VectorWithOpMutable *diag = diag_.get();
	THROW_EXCEPTION(
		!diag, std::logic_error
		,"MatrixSymDiagonalStd::diag(): Error, the diagonal vector has not been set! " );
	return *diag;;
}

const VectorSpace::vec_mut_ptr_t&
MatrixSymDiagonalStd::diag_ptr() const
{
	return diag_;
}

// Overridden from MatrixBase

size_type MatrixSymDiagonalStd::rows() const
{
	return diag_.get() ? diag_->dim() : 0;
}

size_type MatrixSymDiagonalStd::nz() const
{
	return diag_.get() ? diag_->nz() : 0;
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixSymDiagonalStd::space_cols() const {
	return diag_->space();
}

const VectorSpace& MatrixSymDiagonalStd::space_rows() const {
	return diag_->space();
}

MatrixWithOp&
MatrixSymDiagonalStd::operator=(const MatrixWithOp& M)
{
	const MatrixSymDiagonalStd
		*p_M = dynamic_cast<const MatrixSymDiagonalStd*>(&M);

	THROW_EXCEPTION(
		p_M == NULL, std::logic_error
		,"MatrixSymDiagonalStd::operator=(M): Error, the matrix M with concrete type "
		"\'" << typeid(M).name() << "\' does not support the MatrixSymDiagonalStd type! " );

	if( p_M == this ) return *this; // Assignment to self

	diag_    = p_M->diag_;  // lazy copy!
	unique_  = p_M->unique_;

	return *this;
}

bool MatrixSymDiagonalStd::Mp_StM(
	MatrixWithOp* M_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
	// ToDo: validate the vector spaces for the matrices!
	MultiVectorMutable
		*M_mv_lhs = dynamic_cast<MultiVectorMutable*>(M_lhs);
	if(!M_mv_lhs)
		return false;
	VectorSpace::vec_mut_ptr_t
		M_diag = M_mv_lhs->diag(0);
	if(!M_diag.get())
		return false; // Access to the diagonal is not supported!
	Vp_StV( M_diag.get(), alpha, *diag_ );
	return true;
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorWithOp& v_rhs2, value_type beta) const
{
	// ToDo: Validate input!
	if(beta == 0.0)
		*v_lhs = 0.0;
	else if(beta != 1.0)
		Vt_S( v_lhs, beta );
	ele_wise_prod( alpha, v_rhs2, *diag_, v_lhs );
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	MatrixWithOp::Vp_StMtV(v_lhs,alpha,trans_rhs1,sv_rhs2,beta); // ToDo: Implement specialized!
}

// Overridden from MatrixNonsingular

void MatrixSymDiagonalStd::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const VectorWithOp& v_rhs2) const
{
	ele_wise_divide( 1.0, v_rhs2, *diag_, v_lhs );
}

void MatrixSymDiagonalStd::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	MatrixNonsingular::V_InvMtV(v_lhs,trans_rhs1,sv_rhs2 ); // ToDo: Implement specialized!
}

bool MatrixSymDiagonalStd::syrk(
	BLAS_Cpp::Transp   A_trans
	,value_type        a
	,value_type        b
	,MatrixSymWithOp   *B
	) const
{
	MatrixSymDiagonalStd    *B_sd = dynamic_cast<MatrixSymDiagonalStd*>(B);
	if(!B_sd) return false;
	VectorWithOpMutable     &B_diag = B_sd->diag();
	const VectorWithOp      &A_diag = this->diag();
	// B = b*B + a*A*A
	Vt_S( &B_diag, b );
	ele_wise_prod( 1.0, A_diag, A_diag, &B_diag );   // B.diag(i) += a * (A.diag)(i) * (A.diag)(i)
	return true;
}

// Overridden from MatrixSymInitDiagonal

void MatrixSymDiagonalStd::init_identity( const VectorSpace& space_diag, value_type alpha )
{
	diag_ = space_diag.create_member();
	if( diag_->dim() )
		*diag_ = alpha;
}

void MatrixSymDiagonalStd::init_diagonal( const VectorWithOp& diag )
{
	diag_ = diag.space().create_member();
	*diag_ = diag;
}

// Overridden from MatrixSymDiagonal

const VectorWithOp& MatrixSymDiagonalStd::diag() const
{
	return const_cast<MatrixSymDiagonalStd*>(this)->diag();
}

// private

void MatrixSymDiagonalStd::copy_unique()
{
	if( diag_.get() && diag_.count() > 1 && unique_ )
		diag_ = diag_->clone();
}

} // end namespace AbstractLinAlgPack
