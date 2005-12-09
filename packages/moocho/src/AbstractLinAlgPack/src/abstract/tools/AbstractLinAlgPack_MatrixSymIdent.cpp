// ///////////////////////////////////////////////////////////
// MatrixSymIdent.cpp
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

#include <ostream>

#include "AbstractLinAlgPack_MatrixSymIdent.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace AbstractLinAlgPack {

// Constructors/initalizers

MatrixSymIdent::MatrixSymIdent(
	const VectorSpace::space_ptr_t&          vec_space
	,const value_type                        scale
	)
{
	this->initialize(vec_space,scale);
}

void MatrixSymIdent::initialize(
	const VectorSpace::space_ptr_t&          vec_space
	,const value_type                        scale
	)
{
	vec_space_ = vec_space;
	scale_     = scale;
}

// Overridden from MatrixBase

size_type MatrixSymIdent::rows() const
{
	return vec_space_.get() ? vec_space_->dim() : 0;
}

size_type MatrixSymIdent::nz() const
{
	return vec_space_.get() ? vec_space_->dim() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixSymIdent::space_cols() const {
	return *vec_space_;
}

std::ostream& MatrixSymIdent::output(std::ostream& out) const
{
	out << "Identity matrix of dimension " << rows() << " x " << rows() << std::endl;
	return out;
}

void MatrixSymIdent::Vp_StMtV(
	VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	,const Vector& x, value_type b
	) const
{
	AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
	Vt_S(y,b);
    Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsing

void MatrixSymIdent::V_InvMtV(
	VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
	) const
{
	AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
	LinAlgOpPack::V_StV(y,1.0/scale_,x);
}

} // end namespace AbstractLinAlgPack
