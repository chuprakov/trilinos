// //////////////////////////////////////////////////////////////////
// MatrixSymWithOpNonsingularSerial.h
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

#ifndef SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H

#include "MatrixSymWithOpSerial.h"
#include "MatrixSymNonsingularSerial.h"
#include "AbstactLinAlgPack/include/MatrixWithOpNonsingular.h"

namespace SparseLinAlgPack {

///
/** Abstract base class for all serial polymorphic symmetric nonsingular matrices that
  * can be used to compute matrix-vector products and solve for linear systems relatively
  * efficiently.
  */
class MatrixSymWithOpNonsingularSerial 
	: virtual public MatrixSymWithOp
	, virtual public MatrixSymNonsingular
	, virtual public MatrixWithOpNonsingular
	, virtual public MatrixSymWithOpNonsingular
{};

}	// end namespace SparseLinAlgPack

#endif	// SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H
