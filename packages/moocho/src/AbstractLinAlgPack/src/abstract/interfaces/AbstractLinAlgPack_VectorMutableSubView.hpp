// ////////////////////////////////////////////////////////////
// VectorWithOpMutableSubView.hpp
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

#include "VectorWithOpMutable.hpp"
#include "VectorWithOpSubView.hpp"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a sub-view of a VectorWithOpMutable object.
 *
 * Not all of the methods from VectorWithOpMutable are overridden, only those that
 * need to be or may result in better performance.
 *
 * The default constructor and copy constructors are allowd but the default assignment
 * operator is not allowed since it does not have the correct sematics.
 *
 * There is really not much to this vector subclass.  The subclass is only possible
 * because of the \c first_ele, \c sub_dim, and \c global_offset options with apply_transforamtion().  The
 * vector space object returned by <tt>this->space()</tt> is of type \c VectorSpaceSubSpace
 * which in turn relys on \c VectorSpace::sub_space().
 */
class VectorWithOpMutableSubView
	: virtual public VectorWithOpMutable
	, virtual public VectorWithOpSubView
{
public:

	///
	/** Constructs to uninitialized.
	 *
	 * Postconditions: see \c set_uninitialized().
	 */
	VectorWithOpMutableSubView();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorWithOpMutableSubView( const vec_mut_ptr_t& full_vec, const Range1D& rng );

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

	/** @name Overridden from VectorWithOp */
	//@{

	/// Overridden to pick VectorWithOpSubView::sub_view().
	vec_ptr_t sub_view( const Range1D& rng ) const;

	//@}

	/** @name Overridden from VectorWithOpMutable */
	//@{
	
	///
	/** Calls \c apply_transformation() on the underlying full vectors.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>dynamic_cast<const VectorWithOpSubView*>(vecs[k]) != NULL</tt>, for <tt>k=0..num_vecs</tt>
	 *      (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>dynamic_cast<VectorWithOpMutableSubView*>(targ_vecs[k]) != NULL</tt>, for <tt>k=0..num_targ_vecs</tt>
	 *      (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>dynamic_cast<const VectorWithOpSubView*>(vecs[k])->full_vec()->space().is_compatible(
	 *      this->full_vec()->space() ) == true</tt>, for <tt>k=0..num_vecs</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>dynamic_cast<VectorWithOpMutableSubView>(targ_vecs[k])->full_vec()->space().is_compatible(
	 *      this->full_vec()->space() ) == true</tt>, for <tt>k=0..num_targ_vecs</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 */
	void apply_transformation(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		);
	///
	void set_ele( index_type i, value_type val );
	///
	vec_mut_ptr_t sub_view( const Range1D& rng );
	///
	void set_sub_vector( const RTOp_SubVector& sub_vec );

	//@}

private:

	vec_mut_ptr_t       full_vec_; // The full vector

	// Not defined and not to be called!
	VectorWithOpMutableSubView& operator=(const VectorWithOpMutableSubView&);

}; // end class VectorWithOpMutableSubView

// //////////////////////////////////
// Inline members

inline
VectorWithOpMutableSubView::VectorWithOpMutableSubView()
{}

inline
const VectorWithOpMutableSubView::vec_mut_ptr_t&
VectorWithOpMutableSubView::full_vec() const
{
	return full_vec_;
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
