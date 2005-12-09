// //////////////////////////////////////////////////////////
// MatrixSymOp.cpp
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

#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"

namespace AbstractLinAlgPack {

MatrixSymOp::mat_mswo_mut_ptr_t
MatrixSymOp::clone_mswo()
{
	return Teuchos::null;
}

MatrixSymOp::mat_mswo_ptr_t
MatrixSymOp::clone_mswo() const
{
	return Teuchos::null;
}

void MatrixSymOp::Mp_StPtMtP(
	MatrixSymOp* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

void MatrixSymOp::Mp_StMtMtM(
	MatrixSymOp* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

// Overridden from MatrixOp


size_type MatrixSymOp::cols() const
{
	return this->rows();
}

const VectorSpace& MatrixSymOp::space_rows() const
{
	return this->space_cols();
}

MatrixSymOp::mat_mut_ptr_t
MatrixSymOp::clone()
{
	return clone_mswo();
}

MatrixSymOp::mat_ptr_t
MatrixSymOp::clone() const
{
	return clone_mswo();
}

}	// end namespace AbstractLinAlgPack 
