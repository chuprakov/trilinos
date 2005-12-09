// ///////////////////////////////////////
// AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp
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

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

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
	  * Warning, the entire DMatrixSlice InvChol->gms() can be
	  * used for workspace!
	  */
	virtual void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const = 0;	

};	// end class MatrixExtractInvCholFactor

}	// end namespace AbstractLinAlgPack

#endif // MATRIX_EXTRACT_INV_CHOL_FACTOR_H
