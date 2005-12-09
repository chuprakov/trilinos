// ///////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymInitDiag.hpp
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

#ifndef MATRIX_SYM_INIT_DIAGONAL_H
#define MATRIX_SYM_INIT_DIAGONAL_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Mix-in Interface for setting a matrix to a diagonal {abstract}.
 */
class MatrixSymInitDiag {
public:

	///
	virtual ~MatrixSymInitDiag() {}

	/// Initialize a <tt>n x n</tt> identity matrix scaled by \c alpha (where <tt>n = diag.dim()</tt>).
	virtual void init_identity( const VectorSpace& space_diag, value_type alpha = 1.0 ) = 0;

	/// Initialize an <tt>n x n</tt> diagonal matrix (where <tt>n = diag.dim()</tt>).
	virtual void init_diagonal( const Vector& diag ) = 0;

}; // end class MatrixSymInitDiag

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_INIT_DIAGONAL_H
