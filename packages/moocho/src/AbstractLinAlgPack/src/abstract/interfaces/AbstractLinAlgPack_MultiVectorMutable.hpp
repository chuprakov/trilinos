// ///////////////////////////////////////////////////////////
// MultiVectorMutable.h
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

#ifndef MULTI_VECTOR_MUTABLE_H
#define MULTI_VECTOR_MUTABLE_H

#include "MultiVector.h"

namespace AbstractLinAlgPack {

///
/** Mix-in interface of providing mutable row/column/diagonal access to a matrix object.
 *
 * This interface extends the \c MutiVector interface and is ment to be included
 * in a subclass alone with \c MatrixWithOp in order to add this very specialized
 * functionality.
 *
 * These vectors allow the modification of the matrix
 * row by row and column by column.  Each of the views
 * is transient and should be ussed and then discarded.
 *
 * Note that the matrix is only guaranteed to be modified
 * after the smart reference counted pointer returned from
 * these methods is destoryed.  Consider the following code:
 \verbatim

 void f( MultiVectorMutable* M, index_type i )
 {
	MultiVectorMutable::vec_mut_ptr_t
	    row_i =M->row(i);
	*row_i = 0.0;
	// The underlying matrix may not be modified at this point.
	row_i = NULL;
	// Now the underlying matrix is guaranteed to be modified and
	// we can assume this in the following code.
	...
 }
 \endverbatim
 *
 * Default implementations of the const methods from \c MultiVector call the
 * non-const methods defined here and cast the pointers.
 *
 * ToDo: Finish documentation!
 */
class MultiVectorMutable : virtual public MultiVector {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<VectorWithOpMutable>  vec_mut_ptr_t;

	/** @name Provide mutable row, column and/or diagonal access. */
	//@{

	///
	/** Get a mutable row vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_ROW</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t row(index_type i) = 0;
	///
	/** Get a mutable column vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_COL</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t col(index_type j) = 0;
	///
	/** Get a mutable diagonal vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_DIAG</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t diag(int k) = 0;

	//@}

	// ToDo: Add apply_transformation(...) to operate over all of the constituent vectors.
	// Note, the client must also supply a reduction operator object that can reduce
	// reduction object completled across rows or columns etc.  The default implementation
	// can be given here in terms of ros(), col(), or diag().

	/** @name Overridden from MultiVector */
	//@{

	///
	virtual vec_ptr_t row(index_type i) const;
	///
	virtual vec_ptr_t col(index_type j) const;
	///
	virtual vec_ptr_t diag(int k) const;

	//@}

}; // end class MultiVectorMutable

} // end namespace AbstractLinAlgPack

#endif // MULTI_VECTOR_MUTABLE_H
