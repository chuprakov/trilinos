// //////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingular.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H

#include "MatrixWithOp.h"
#include "MatrixNonsingular.h"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all nonsingular polymorphic matrices
 * that can be used to compute matrix-vector products and solve
 * for linear systems efficiently.
 */
class MatrixWithOpNonsingular
	: virtual public MatrixWithOp
	, virtual public MatrixNonsingular
{};

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H
