// ////////////////////////////////////////////////////////////
// VectorWithOpMutableSubView.h
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

#include "VectorWithOpMutable.h"
#include "VectorWithOpSubView.h"

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

	/** @name Overridden from VectorWithOp */
	//@{

	/// Overridden to pick VectorWithOpSubView::sub_view().
	vec_ptr_t sub_view( const Range1D& rng ) const;

	//@}

	/** @name Overridden from VectorWithOpMutable */
	//@{
	
	///
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

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
