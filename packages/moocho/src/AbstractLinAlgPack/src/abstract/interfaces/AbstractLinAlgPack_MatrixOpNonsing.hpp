// //////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingular.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H

#include "MatrixWithOp.h"
#include "MatrixNonsingular.h"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all nonsingular polymorphic matrices
 * that can be used to compute matrix-vector products and solve
 * for linear systems efficiently.
 */
class MatrixWithOpNonsingular
	: virtual public MatrixWithOp
	, virtual public MatrixNonsingular
{
public:

	/** @name Public types */
	//@{

#ifndef DOXYGEN_COMPILE
	///
	typedef MemMngPack::ref_count_ptr<const MatrixWithOpNonsingular>    mat_mwons_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<MatrixWithOpNonsingular>          mat_mwons_mut_ptr_t;
#endif

	//@}

	/** @name Clone */
	//@{

	///
	/** Clone the non-const matrix object (if supported).
	 *
	 * The default implementation returns NULL which is perfectly acceptable.
	 * A matrix object is not required to return a non-NULL value but almost
	 * every good matrix implementation will.
	 */
	virtual mat_mwons_mut_ptr_t clone_mwons();

	///
	/** Clone the const matrix object (if supported).
	 *
	 * The behavior of this method is the same as for the non-const version
	 * above except it returns a smart pointer to a const matrix object.
	 */
	virtual mat_mwons_ptr_t clone_mwons() const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{
	/// Returns <tt>this->clone_mwons()</tt>.
	mat_mut_ptr_t clone();
	/// Returns <tt>this->clone_mwons()</tt>.
	mat_ptr_t clone() const;
	//@}

	/** @name Overridden from MatrixNonsingular */
	//@{
	/// Returns <tt>this->clone_mwons()</tt>.
	mat_mns_mut_ptr_t clone_mns();
	/// Returns <tt>this->clone_mwons()</tt>.
	mat_mns_ptr_t clone_mns() const;
	//@}

	/// Calls operator=(MatrixWithOp&)
	MatrixWithOpNonsingular& operator=(const MatrixWithOpNonsingular& M)
	{ static_cast<MatrixWithOp*>(this)->operator=(M); return *this; }

}; // end class MatrixWithOpNonsingular

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_NONSINGULAR_H
