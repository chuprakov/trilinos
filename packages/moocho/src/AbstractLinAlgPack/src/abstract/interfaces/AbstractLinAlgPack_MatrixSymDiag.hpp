// ///////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymDiag.hpp
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

#ifndef MATRIX_SYM_DIAGONAL_H
#define MATRIX_SYM_DIAGONAL_H

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace AbstractLinAlgPack {

///
/** Interface to all diagonal matrices {abstract}.
 */
class MatrixSymDiag
	: public virtual MatrixSymOpNonsing
{
public:

	/// Give const access to the diagonal
	virtual const Vector& diag() const = 0;

}; // end class MatrixSymDiag

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_DIAGONAL_H
