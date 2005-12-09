// /////////////////////////////////////////////////////////////////////////////
// MatrixSymNonsing.cpp
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

#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"

namespace AbstractLinAlgPack {

MatrixSymNonsing::mat_msns_mut_ptr_t
MatrixSymNonsing::clone_msns()
{
	return Teuchos::null;
}

MatrixSymNonsing::mat_msns_ptr_t
MatrixSymNonsing::clone_msns() const
{
	return Teuchos::null;
}

void MatrixSymNonsing::M_StMtInvMtM(
	  MatrixSymOp* S, value_type a, const MatrixOp& B
	, BLAS_Cpp::Transp B_trans, EMatrixDummyArg ) const
{
	assert(0); // ToDo: Implement!
}

// Overridden from MatrixNonsing

MatrixSymNonsing::mat_mns_mut_ptr_t
MatrixSymNonsing::clone_mns()
{
	return clone_msns();
}

MatrixSymNonsing::mat_mns_ptr_t
MatrixSymNonsing::clone_mns() const
{
	return clone_msns();
}

}	// end namespace AbstractLinAlgPack
