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

#include "AbstractLinAlgPack/include/MatrixSymWithOp.h"
#include "AbstractLinAlgPack/include/EtaVector.h"

namespace AbstractLinAlgPack {

size_type MatrixSymWithOp::cols() const
{
	return this->rows();
}

void MatrixSymWithOp::Mp_StPtMtP(
	MatrixSymWithOpMutable* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

void MatrixSymWithOp::Mp_StMtMtM(
	MatrixSymWithOpMutable* sym_lhs, value_type alpha
	, EMatRhsPlaceHolder dummy_place_holder
	, const MatrixWithOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

}	// end namespace AbstractLinAlgPack 
