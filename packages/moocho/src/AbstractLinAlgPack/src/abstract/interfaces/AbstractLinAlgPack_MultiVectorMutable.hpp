// /////////////////////////////////////////////////////////////////////////////
// MultiVectorMutable.hpp
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

#ifndef ALAP_MULTI_VECTOR_MUTABLE_H
#define ALAP_MULTI_VECTOR_MUTABLE_H

#include "MultiVector.hpp"

namespace AbstractLinAlgPack {

///
/** Interface for a collection of mutable vectors (multi-vector, matrix).
 *
 * This interface extends the \c MutiVector interface an allows mutable access to
 * the constituent vectors.
 *
 * These vectors allow the modification of the matrix row by row, column by column,
 * and/or diagonal by diagonal.  Each of the views is transient and should be used
 * and discarded quickly.
 *
 * Note that the underlying matrix is only guaranteed to be modified after the smart reference
 * counted pointer returned from these methods is destoryed.  For example, consider the following code:
 \code

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
 \endcode
 * Default implementations of the const access methods \c row() \c col()
 * and \c diag() from \c MultiVector call the non-const methods defined
 * here and cast the pointers.
 *
 * Many of the default implementations of the linear algebra operations in
 * \c MatrixOp and the other matrix interfaces rely on the left hand side
 * matrix objects supporting the \c MultiVectorMutable interface.
 */
class MultiVectorMutable : virtual public MultiVector {
public:
	
	///
	using MultiVector::col;
	///
	using MultiVector::row;
	///
	using MultiVector::diag;

	///
	typedef MemMngPack::ref_count_ptr<VectorMutable>         vec_mut_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MultiVectorMutable>    multi_vec_mut_ptr_t;

	/** @name Provide mutable row, column and/or diagonal access */
	//@{

	///
	/** Get a mutable column vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & COL_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>return.get() != NULL</tt>] <tt>space_cols().is_compatible(return->space()) == true</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t col(index_type j) = 0;
	///
	/** Get a mutable row vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & ROW_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>return.get() != NULL</tt>] <tt>space_rows().is_compatible(return->space()) == true</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t row(index_type i) = 0;
	///
	/** Get a mutable diagonal vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & DIAG_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual vec_mut_ptr_t diag(int k) = 0;

	//@}

	/** @name Sub-view methods */
	//@{

	///
	/** Returns a mutable sub-view of the multi vector.
	 *
	 * ToDo: Finish documentation!
	 *
	 * The default implementation returns a \c MultiVectorMutableSubView object for
	 * any valid arbitary sub-view.
	 */
	virtual multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);
	
	///
	/** Inlined implementation calls <tt>this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu))</tt>.
	 */
	multi_vec_mut_ptr_t mv_sub_view(
		const index_type& rl, const index_type& ru
		,const index_type& cl, const index_type& cu
		);

	//@}

	/** @name Overridden from MultiVector */
	//@{

	///
	virtual vec_ptr_t col(index_type j) const;
	///
	virtual vec_ptr_t row(index_type i) const;
	///
	virtual vec_ptr_t diag(int k) const;
	///
	multi_vec_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const;

	//@}

}; // end class MultiVectorMutable

// //////////////////////////////////////////////////
// Inlined methods for MultiVector

inline
MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(
	const index_type& rl, const index_type& ru
	,const index_type& cl, const index_type& cu
	)
{
	return this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu));
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_MUTABLE_H
