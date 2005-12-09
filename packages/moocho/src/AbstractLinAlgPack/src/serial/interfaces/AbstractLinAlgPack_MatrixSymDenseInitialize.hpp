// //////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymDenseInitialize.hpp
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

#ifndef MATRIX_SYM_DENSE_INITIALIZE_H
#define MATRIX_SYM_DENSE_INITIALIZE_H

#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

///
/** Mix-in Interface for initializing a matrix with a dense symmetric matrix.
  */
class MatrixSymDenseInitialize
	: virtual public AbstractLinAlgPack::MatrixBase // doxygen needs the full name
{
public:

	///
	/** Initialize with a symmetric dense matrix.
	 *
	 * Through this interface there are absolutly no postconditions
	 * as the the state of \c this after this function executes.
	 * The implementation can use \c M to initialize itself any way
	 * it would like.
	 */
	virtual void initialize( const DMatrixSliceSym& M ) = 0;

};	// end class MatrixSymDenseInitialize

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_SYM_DENSE_INITIALIZE_H
