// /////////////////////////////////////////////////////////////////////////////
// MatrixSymNonsingular.cpp
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

#include "AbstractLinAlgPack/src/MatrixSymNonsingular.h"
#include "AbstractLinAlgPack/src/EtaVector.h"

namespace AbstractLinAlgPack {

MatrixSymNonsingular::mat_msns_mut_ptr_t
MatrixSymNonsingular::clone_msns()
{
	return MemMngPack::null;
}

MatrixSymNonsingular::mat_msns_ptr_t
MatrixSymNonsingular::clone_msns() const
{
	return MemMngPack::null;
}

void MatrixSymNonsingular::M_StMtInvMtM(
	  MatrixSymWithOp* S, value_type a, const MatrixWithOp& B
	, BLAS_Cpp::Transp B_trans, EMatrixDummyArg ) const
{
	assert(0); // ToDo: Implement!
}

// Overridden from MatrixNonsingular

MatrixSymNonsingular::mat_mns_mut_ptr_t
MatrixSymNonsingular::clone_mns()
{
	return clone_msns();
}

MatrixSymNonsingular::mat_mns_ptr_t
MatrixSymNonsingular::clone_mns() const
{
	return clone_msns();
}

}	// end namespace AbstractLinAlgPack
