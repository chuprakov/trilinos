// //////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingularSerial.hpp
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

#ifndef SLAP_MATRIX_WITH_OP_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_WITH_OP_NONSINGULAR_SERIAL_H

#include "MatrixOp.hpp"
#include "AbstractLinAlgPack/src/MatrixNonsing.hpp"

namespace SparseLinAlgPack {

///
/** Abstract base class for all serial nonsingular polymorphic matrices
 * that can be used to compute matrix-vector products and solve for
 * linear systems efficiently.
 */
class MatrixWithOpNonsingularSerial
	: virtual public MatrixWithOpSerial
	, virtual public MatrixNonsingularSerial
	, virtual public AbstractLinAlgPack::MatrixOpNonsing  // doxygen needs full name
{};

}	// end namespace SparseLinAlgPack

#endif	// SLAP_MATRIX_WITH_OP_NONSINGULAR_SERIAL_H
