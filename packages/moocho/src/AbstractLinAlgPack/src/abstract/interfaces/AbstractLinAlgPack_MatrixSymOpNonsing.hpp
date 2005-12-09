// //////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymOpNonsing.hpp
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H

#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all polymorphic symmetrix nonsingular matrices that
 * can be used to compute matrix-vector products and solve for
 * linear systems relatively efficently.
 */
class MatrixSymOpNonsing 
	: virtual public MatrixSymOp
	, virtual public MatrixSymNonsing
	, virtual public MatrixOpNonsing
{
public:

	/** @name Public types */
	//@{

#ifndef DOXYGEN_COMPILE
	///
	typedef Teuchos::RefCountPtr<const MatrixSymOpNonsing>    mat_mswons_ptr_t;
	///
	typedef Teuchos::RefCountPtr<MatrixSymOpNonsing>          mat_mswons_mut_ptr_t;
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
	virtual mat_mswons_mut_ptr_t clone_mswons();

	///
	/** Clone the const matrix object (if supported).
	 *
	 * The behavior of this method is the same as for the non-const version
	 * above except it returns a smart pointer to a const matrix object.
	 */
	virtual mat_mswons_ptr_t clone_mswons() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mut_ptr_t clone();
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_ptr_t clone() const;
	//@}

	/** @name Overridden from MatrixNonsing */
	//@{
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mns_mut_ptr_t clone_mns();
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mns_ptr_t clone_mns() const;
	//@}

	/** @name Overridden from MatrixSymOp */
	//@{
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mswo_mut_ptr_t clone_mswo();
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mswo_ptr_t clone_mswo() const;
	//@}

	/** @name Overridden from MatrixSymNonsing */
	//@{
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_msns_mut_ptr_t clone_msns();
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_msns_ptr_t clone_msns() const;
	//@}

	/** @name Overridden from MatrixOpNonsing */
	//@{
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mwons_mut_ptr_t clone_mwons();
	/// Returns <tt>this->clone_mswons()</tt>.
	mat_mwons_ptr_t clone_mwons() const;
	//@}

	/// Calls operator=(MatrixOp&)
	MatrixSymOpNonsing& operator=(const MatrixSymOpNonsing& M)
	{ static_cast<MatrixOp*>(this)->operator=(M); return *this; }

}; // end class MatrixSymOpNonsing

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H
