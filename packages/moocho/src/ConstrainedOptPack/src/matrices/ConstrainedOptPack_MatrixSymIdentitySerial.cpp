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

#include <limits>

#include "ConstrainedOptPack_MatrixSymIdentitySerial.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace ConstrainedOptPack {

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

// Overridden from MatrixOp

std::ostream& MatrixSymIdentitySerial::output(std::ostream& out) const
{
	out << "Identity matrix of size " << size_ << " x " << size_ << std::endl;
	return out;
}

// Overridden from MatrixOpSerial

void MatrixSymIdentitySerial::Vp_StMtV(
	DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	,const DVectorSlice& x, value_type b
	) const
{
	DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
	DenseLinAlgPack::Vt_S(y,b);
	DenseLinAlgPack::Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsinguarSerial

void MatrixSymIdentitySerial::V_InvMtV(
	DVectorSlice* y, BLAS_Cpp::Transp M_trans, const DVectorSlice& x
	) const
{
	DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
	LinAlgOpPack::V_StV(y,scale_,x);
}

// Overridden from MatrixSymNonsing

void MatrixSymIdentitySerial::M_StMtInvMtM(
	  DMatrixSliceSym* S, value_type a
	  ,const MatrixOpSerial& B, BLAS_Cpp::Transp B_trans
	  ,EMatrixDummyArg dummy_arg
	) const
{
	this->MatrixSymNonsingSerial::M_StMtInvMtM(S,a,B,B_trans,dummy_arg);
	// ToDo: Implement by calling S = b*S + scale*a*op(B')*op(B)
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymIdentitySerial::extract_inv_chol( DMatrixSliceTriEle* InvChol ) const
{
	if( scale_ < 0.0 )
		throw std::logic_error(
			"MatrixSymIdentitySerial::extract_inv_chol(...) : "
			"Error, we can not compute the inverse cholesky factor "
			"of a negative definite matrix." );
	DenseLinAlgPack::assign( &InvChol->gms(), 0.0 );
	InvChol->gms().diag() = 1.0 / ::sqrt(scale_);
}

} // end namespace ConstrainedOptPack
