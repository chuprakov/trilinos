// //////////////////////////////////////////////////////////
// MatrixSymWithOp.cpp
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

#include "AbstractLinAlgPack/src/MatrixSymWithOp.hpp"
#include "AbstractLinAlgPack/src/EtaVector.hpp"

namespace AbstractLinAlgPack {

MatrixSymWithOp::mat_mswo_mut_ptr_t
MatrixSymWithOp::clone_mswo()
{
	return MemMngPack::null;
}

MatrixSymWithOp::mat_mswo_ptr_t
MatrixSymWithOp::clone_mswo() const
{
	return MemMngPack::null;
}

void MatrixSymWithOp::Mp_StPtMtP(
	MatrixSymWithOp* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

void MatrixSymWithOp::Mp_StMtMtM(
	MatrixSymWithOp* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const MatrixWithOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

// Overridden from MatrixWithOp


size_type MatrixSymWithOp::cols() const
{
	return this->rows();
}

const VectorSpace& MatrixSymWithOp::space_rows() const
{
	return this->space_cols();
}

MatrixSymWithOp::mat_mut_ptr_t
MatrixSymWithOp::clone()
{
	return clone_mswo();
}

MatrixSymWithOp::mat_ptr_t
MatrixSymWithOp::clone() const
{
	return clone_mswo();
}

}	// end namespace AbstractLinAlgPack 
