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
#include "AbstractLinAlgPack/include/MatrixWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"

namespace AbstractLinAlgPack {

MatrixSymDiagonalStd::MatrixSymDiagonalStd()
{}

VectorWithOpMutable& MatrixSymDiagonalStd::diag()
{
	return *diag_;
}

const VectorWithOp& MatrixSymDiagonalStd::diag() const
{
	return *diag_;
}

// Overridden from Matrix

size_type MatrixSymDiagonalStd::rows() const
{
	return diag_->dim();
}

size_type MatrixSymDiagonalStd::nz() const
{
	return diag_->nz();
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixSymDiagonalStd::space_rows() const {
	return diag_->space();
}

const VectorSpace& MatrixSymDiagonalStd::space_cols() const {
	return diag_->space();
}

void MatrixSymDiagonalStd::Mp_StM(
	MatrixWithOp* M_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
	MatrixWithOp::Mp_StM(M_lhs,alpha,trans_rhs); // ToDo: Implement specialized!
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorWithOp& v_rhs2, value_type beta) const
{
	ele_wise_prod( alpha, v_rhs2, *diag_, v_lhs );
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	MatrixWithOp::Vp_StMtV(v_lhs,alpha,trans_rhs1,sv_rhs2,beta); // ToDo: Implement specialized!
}

// Overridden from MatrixFactorized

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
	MatrixFactorized::V_InvMtV(v_lhs,trans_rhs1,sv_rhs2 ); // ToDo: Implement specialized!
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
