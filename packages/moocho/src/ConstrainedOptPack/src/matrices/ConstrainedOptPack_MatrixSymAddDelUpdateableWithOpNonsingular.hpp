// ////////////////////////////////////////////////////
// ConstrainedOptPack_MatrixSymAddDelUpdateableWithOpNonsingular.hpp
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

#ifndef MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H
#define MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Interface for updating a symmetric matrix and its factorization 
 * by adding and deleting rows and columns and preforming operations with it.
 *
 * ToDo: Finish documentation.
 */
class MatrixSymAddDelUpdateableWithOpNonsingular {
public:

	///
	virtual ~MatrixSymAddDelUpdateableWithOpNonsingular() {}
	///
	virtual const MatrixSymOpNonsing& op_interface() const = 0;
	///
	virtual MatrixSymAddDelUpdateable& update_interface() = 0;
	///
	virtual const MatrixSymAddDelUpdateable& update_interface() const = 0;

};	// end class MatrixSymAddDelUpdateableWithOpNonsingular

}	// namespace ConstrainedOptPack 

#endif	// MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H
