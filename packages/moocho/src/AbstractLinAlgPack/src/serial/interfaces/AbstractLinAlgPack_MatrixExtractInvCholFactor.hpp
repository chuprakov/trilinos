// ///////////////////////////////////////
// MatrixExtractInvCholFactor.hpp
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

#ifndef MATRIX_EXTRACT_INV_CHOL_FACTOR_H
#define MATRIX_EXTRACT_INV_CHOL_FACTOR_H

#include "SparseLinAlgPackTypes.hpp"
#include "AbstractLinAlgPack/src/MatrixBase.hpp"

namespace SparseLinAlgPack {

///
/** Mix-in Interface for extracting the inverse cholesky factor of a dense symmetric
  * positive definite matrix.
  */
class MatrixExtractInvCholFactor
	: virtual public AbstractLinAlgPack::MatrixBase // doxygen needs full name
{
public:

	///
	/** Extract the inverse cholesly factor.
	  *
	  * Warning, the entire GenMatrixSlice InvChol->gms() can be
	  * used for workspace!
	  */
	virtual void extract_inv_chol( tri_ele_gms* InvChol ) const = 0;	

};	// end class MatrixExtractInvCholFactor

}	// end namespace SparseLinAlgPack

#endif // MATRIX_EXTRACT_INV_CHOL_FACTOR_H
