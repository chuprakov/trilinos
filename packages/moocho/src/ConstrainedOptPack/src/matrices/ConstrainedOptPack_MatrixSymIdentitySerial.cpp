// ///////////////////////////////////////////////////////////
// MatrixSymIdentitySerial.cpp
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

#include "ConstrainedOptimizationPack/src/MatrixSymIdentitySerial.hpp"
#include "LinAlgPack/src/GenMatrixAsTriSym.hpp"
#include "LinAlgPack/src/GenMatrixOp.hpp"
#include "LinAlgPack/src/GenMatrixOut.hpp"
#include "LinAlgPack/src/LinAlgOpPack.hpp"
#include "LinAlgPack/src/LinAlgPackAssertOp.hpp"

namespace ConstrainedOptimizationPack {

// Constructors

MatrixSymIdentitySerial::MatrixSymIdentitySerial(size_type size, value_type scale)
{
	this->initialize(size,scale);
}

void MatrixSymIdentitySerial::initialize(size_type size, value_type scale)
{
	size_  = size;
	scale_ = scale;
}

// Overridden from MatrixBase

size_type MatrixSymIdentitySerial::rows() const
{
	return size_;
}

size_type MatrixSymIdentitySerial::nz() const
{
	return size_;
}

// Overridden from MatrixWithOp

std::ostream& MatrixSymIdentitySerial::output(std::ostream& out) const
{
	out << "Identity matrix of size " << size_ << " x " << size_ << std::endl;
	return out;
}

// Overridden from MatrixWithOpSerial

void MatrixSymIdentitySerial::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	,const VectorSlice& x, value_type b
	) const
{
	LinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
	LinAlgPack::Vt_S(y,b);
	LinAlgPack::Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsinguarSerial

void MatrixSymIdentitySerial::V_InvMtV(
	VectorSlice* y, BLAS_Cpp::Transp M_trans, const VectorSlice& x
	) const
{
	LinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
	LinAlgOpPack::V_StV(y,scale_,x);
}

// Overridden from MatrixSymNonsingular

void MatrixSymIdentitySerial::M_StMtInvMtM(
	  sym_gms* S, value_type a
	  ,const MatrixWithOpSerial& B, BLAS_Cpp::Transp B_trans
	  ,EMatrixDummyArg dummy_arg
	) const
{
	this->MatrixSymNonsingularSerial::M_StMtInvMtM(S,a,B,B_trans,dummy_arg);
	// ToDo: Implement by calling S = b*S + scale*a*op(B')*op(B)
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymIdentitySerial::extract_inv_chol( tri_ele_gms* InvChol ) const
{
	if( scale_ < 0.0 )
		throw std::logic_error(
			"MatrixSymIdentitySerial::extract_inv_chol(...) : "
			"Error, we can not compute the inverse cholesky factor "
			"of a negative definite matrix." );
	LinAlgPack::assign( &InvChol->gms(), 0.0 );
	InvChol->gms().diag() = 1.0 / ::sqrt(scale_);
}

} // end namespace ConstrainedOptimizationPack
