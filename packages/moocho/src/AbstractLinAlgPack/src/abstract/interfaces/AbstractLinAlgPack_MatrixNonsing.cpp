// //////////////////////////////////////////////////////////////////////////////////
//  MatrixNonsingular.cpp
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

// ToDo: 3/6/00: Provide default implementations for these
// operations.

#include <assert.h>

#include "AbstractLinAlgPack/include/MatrixNonsingular.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"

namespace AbstractLinAlgPack {

// Level-2 BLAS

void MatrixNonsingular::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	assert(0); // ToDo: Implement!
}

value_type MatrixNonsingular::transVtInvMtV(
	const VectorWithOp& v_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3) const
{
	assert(0); // ToDo: Implement!
}

value_type MatrixNonsingular::transVtInvMtV(
	const SpVectorSlice& sv_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
	assert(0); // ToDo: Implement!
}

// Level-3 BLAS

void MatrixNonsingular::M_StInvMtM(
	MatrixWithOp* m_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
	,BLAS_Cpp::Transp trans_rhs2) const
{
	assert(0); // ToDo: Implement!
}

void MatrixNonsingular::M_StInvMtM(
	MatrixWithOp* g_lhs, value_type alpha
	, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2) const
{
	assert(0); // ToDo: Implement!
}

}	// end namespace AbstractLinAlgPack
