// //////////////////////////////////////////////////////////////////////////////////
// MatrixBase.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_BASE_H
#define ABSTRACT_LINALG_PACK_MATRIX_BASE_H

#include <stdexcept>

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

///
/** Base class for all polymorphic matrices.
  */
class MatrixBase {
public:

	/// Thrown if matrices are incompatible
	class IncompatibleMatrices : public std::logic_error
	{public: IncompatibleMatrices(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Virtual destructor
	virtual ~MatrixBase() {}

	/// Return the number of rows in the matrix
	virtual size_type rows() const = 0;

	/// Return the number of columns in the matrix
	virtual size_type cols() const = 0;

	///
	/** Return the number of nonzero elements in the matrix.
	  *
	  * The default is to just assume it is dense and to return
	  * rows() * cols().
	  */
	virtual size_type nz() const
	{
		return rows() * cols();
	}

};	// end class MatrixBase

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_BASE_H
