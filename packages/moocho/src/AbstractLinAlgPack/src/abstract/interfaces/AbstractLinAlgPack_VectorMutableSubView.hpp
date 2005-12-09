// ////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorMutableSubView.hpp
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

#ifndef VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
#define VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H

#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSubView.hpp"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a sub-view of a VectorMutable object.
 *
 * Not all of the methods from VectorMutable are overridden, only those that
 * need to be or may result in better performance.
 *
 * The default constructor and copy constructors are allowd but the default assignment
 * operator is not allowed since it does not have the correct sematics.
 *
 * There is really not much to this vector subclass.  The subclass is only possible
 * because of the \c first_ele, \c sub_dim, and \c global_offset options with apply_op().  The
 * vector space object returned by <tt>this->space()</tt> is of type \c VectorSpaceSubSpace
 * which in turn relys on \c VectorSpace::sub_space().
 */
class VectorMutableSubView
	: virtual public VectorMutable
	, virtual public VectorSubView
{
public:

	///
	/** Constructs to uninitialized.
	 *
	 * Postconditions: see \c set_uninitialized().
	 */
	VectorMutableSubView();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorMutableSubView( const vec_mut_ptr_t& full_vec, const Range1D& rng );

	///
	/** Initialize.
	 *
	 * Constructs a view of the vector this = vec(rng).
	 *
	 * @param  full_vec  [in] The original full vector.  It is allowed for <tt>full_vec.get() == NULL</tt>
	 *                   in which case <tt>this</tt> is uninitialized (i.e. <tt>this->dim() == 0</tt>).
	 * @param  rng       [in] The range of elements in <tt>full_vec</tt> that <tt>this</tt> vector will represent.
	 */
	void initialize( const vec_mut_ptr_t& vec, const Range1D& rng );

	///
	/** Set uninitialized()
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->dim() == 0</tt>
	 * <li> <tt>this->full_vec() = NULL</tt>
	 * </ul>
	 */
	void set_uninitialized();

	///
	const vec_mut_ptr_t& full_vec() const;

	/** @name Overridden from Vector */
	//@{

	/// Overridden to pick VectorSubView::sub_view().
	vec_ptr_t sub_view( const Range1D& rng ) const;

	//@}

	/** @name Overridden from VectorMutable */
	//@{
	
	///
	void set_ele( index_type i, value_type val );
	///
	vec_mut_ptr_t sub_view( const Range1D& rng );
	///
	void get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec );
	///
	void commit_sub_vector( RTOpPack::MutableSubVector* sub_vec );
	///
	void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );

	//@}

private:

	vec_mut_ptr_t       full_vec_; // The full vector

	// Not defined and not to be called!
	VectorMutableSubView& operator=(const VectorMutableSubView&);

}; // end class VectorMutableSubView

// //////////////////////////////////
// Inline members

inline
VectorMutableSubView::VectorMutableSubView()
{}

inline
const VectorMutableSubView::vec_mut_ptr_t&
VectorMutableSubView::full_vec() const
{
	return full_vec_;
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
