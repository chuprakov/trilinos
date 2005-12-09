// //////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixBase.hpp
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

#include "AbstractLinAlgPack_Types.hpp"

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

	/** @name Vector spaces for the columns and rows of the matrix */
	//@{

	/// Vector space for vectors that are compatible with the columns of the matrix.
	virtual const VectorSpace& space_cols() const = 0;

	/// Vector space for vectors that are compatible with the rows of the matrix.
	virtual const VectorSpace& space_rows() const = 0;

	//@}

	/** @name Dimensionality */
	//@{

	///
	/** Return the number of rows in the matrix.
	 *
	 * The default implementation returns <tt>space_cols().dim()</tt>.
	 */
	virtual size_type rows() const;

	///
	/** Return the number of columns in the matrix.
	 *
	 * The default implementation returns <tt>space_rows().dim()</tt>.
	 */
	virtual size_type cols() const;

	///
	/** Return the number of nonzero elements in the matrix.
	 *
	 * The default is to just assume it is dense and to return
	 * <tt>rows() * cols()</tt>.
	 */
	virtual size_type nz() const;

	//@}

};	// end class MatrixBase

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_BASE_H
