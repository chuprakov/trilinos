// ///////////////////////////////////////////////////////////
// MultiVector.h
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

#ifndef ALAP_MULTI_VECTOR_H
#define ALAP_MULTI_VECTOR_H

#include "MatrixWithOp.h"
#include "ref_count_ptr.h"

namespace AbstractLinAlgPack {

///
/** Matrix interface for a matrix stored as a set of vectors.
 *
 * This interface is quite restrictive in that it allows a client
 * to access a matrix by accessing rows, columns and/or diagonals.
 * The vector objects returned from these access methods are
 * abstract vectors so there is still good implementation flexibility
 * but many matrix implementations will not be able to support
 * this interface.  This is somewhat of a "last resort" interface
 * that allows many matrix operations to have default implementations
 * based on vector operations.
 *
 * Note that only certain kinds of access may be preferred and it is allowed
 * for subclasses to return \c NULL vector objects for some types of access.  For
 * example, a matrix may be naturally oriented by column but row or diagonal
 * access may be very inefficient.  For this reason, the client should call the
 * \c access_by() method which returns a bit field that the client can compare
 * to the constants \c BY_ROW, \c BY_COL and \c BY_DIAG.  The method \c access_by()
 * only returns the types of access that are guarrentted to be efficient, but
 * does not necessarily imply that a type of access is not supported.  For example,
 * <tt>(this->access_by() & BY_ROW) == false</tt> this does not mean that 
 * <tt>this->row(i)</tt> will return \c NULL, it only means that row access will
 * be inefficient.  To determine if a certain type of access is even possible, check
 * the return for \c row(), \c col() and/or \c diag().  For example, if <tt>this->rows(1)</tt>
 * returns \c NULL, then this is a flag that row access is not supported.
 *
 * Note that since, this interface is derived from \c MatrixWithOp that it must
 * support the methods \c space_rows() and \c space_cols().  This does not imply
 * however that the methods \c row() or \c col() must return non-<tt>NULL</tt>.
 *
 * Examples of matrix implementations that can support this interface are a dense
 * BLAS compatible matrix (\c BY_ROW, \c BY_COL and \c BY_DIAG), a compressed column
 * sparse matrix (\c BY_COL only), a compressed row sparse matrix (\c BY_ROW only)
 * etc.
 *
 * ToDo: Finish documentation!
 */
class MultiVector : virtual public MatrixWithOp {
public:

	///
	typedef int  access_by_t;

	///
	enum {
		BY_ROW    = 0x1 ///< 
		,BY_COL   = 0x2 ///<
		,BY_DIAG  = 0x4 ///<
	};

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorWithOp>   vec_ptr_t;

	/** @name Provide row, column and diagonal access as non-mutable vectors */
	//@{

	///
	/** Return a bit field for the types of access that are the most convenient.
	 */
	virtual access_by_t access_by() const = 0;

	///
	/** Get a non-mutable row vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_ROW</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_ptr_t row(index_type i) const = 0;
	///
	/** Get a non-mutable column vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_COL</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_ptr_t col(index_type j) const = 0;
	///
	/** Get a non-mutable diagonal vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & BY_DIAG</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_ptr_t diag(int k) const = 0;

	//@}

	// ToDo: Add apply_reduction(...) to operate over all of the constituent vectors.
	// Note, the client must also supply a reduction operator object that can reduce
	// reduction object completled across rows or columns etc.  The default implementation
	// can be given here in terms of ros(), col(), or diag().

}; // end class MultiVector

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_H
