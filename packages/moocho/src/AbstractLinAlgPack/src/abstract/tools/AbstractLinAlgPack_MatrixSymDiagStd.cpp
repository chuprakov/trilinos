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

#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixSpaceStd.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

MatrixSymDiagonalStd::MatrixSymDiagonalStd( const VectorSpace::vec_mut_ptr_t& diag )
{
	this->initialize(diag);
}

void MatrixSymDiagonalStd::initialize( const VectorSpace::vec_mut_ptr_t& diag )
{
	diag_ = diag;
}

VectorWithOpMutable& MatrixSymDiagonalStd::diag()
{
	VectorWithOpMutable *diag = diag_.get();
	THROW_EXCEPTION(
		!diag, std::logic_error
		,"MatrixSymDiagonalStd::diag(): Error, the diagonal vector has not been set! " );
	return *diag;;
}

const VectorWithOp& MatrixSymDiagonalStd::diag() const
{
	return const_cast<MatrixSymDiagonalStd*>(this)->diag();
}

// Overridden from Matrix

size_type MatrixSymDiagonalStd::rows() const
{
	return diag_.get() ? diag_->dim() : 0;
}

size_type MatrixSymDiagonalStd::nz() const
{
	return diag_.get() ? diag_->nz() : 0;
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixSymDiagonalStd::space_rows() const {
	return diag_->space();
}

const VectorSpace& MatrixSymDiagonalStd::space_cols() const {
	return diag_->space();
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
	assert(0); // ToDo: Finish!
}

void MatrixSymDiagonalStd::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	MatrixNonsingular::V_InvMtV(v_lhs,trans_rhs1,sv_rhs2 ); // ToDo: Implement specialized!
}

// Overridden from MatrixSymInitDiagonal

void MatrixSymDiagonalStd::init_identity( const VectorSpace& space_diag, value_type alpha )
{
	// ToDo: Initalize mat_space_ to a proper matrix space object (need a lot of 
	// new code to create a compatible mutable matrix object).
	diag_ = space_diag.create_member();
	if( diag_->dim() )
		*diag_ = alpha;
}

void MatrixSymDiagonalStd::init_diagonal( const VectorWithOp& diag )
{
	// ToDo: Initalize mat_space_ to a proper matrix space object (need a lot of 
	// new code to create a compatible mutable matrix object).
	diag_  = diag.space().create_member();
	*diag_ = diag;
}

} // end namespace AbstractLinAlgPack
