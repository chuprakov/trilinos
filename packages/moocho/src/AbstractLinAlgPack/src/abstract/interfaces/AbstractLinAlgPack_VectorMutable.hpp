// //////////////////////////////////////////
// VectorWithOpMutable.h
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

#ifndef VECTOR_WITH_OP_MUTABLE_H
#define VECTOR_WITH_OP_MUTABLE_H

#include "VectorBaseMutable.h"
#include "VectorWithOp.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

///
/** \brief Abstract interface for mutable coordinate vectors {abstract}.
  *
  * Objects of this type can act as a target vector of a transformation operation.
  * Similarly to \c VectorWithOp this interface contains very few (only one extra) pure
  * virtual methods that must be overridden.  However, more efficient and more general
  * implementations will choose to override more methods.
  * 
  * The most important method is \c apply_transformation() that allows clients to apply
  * user defined reduction/transformation operators.  Every standard (i.e. BLAS) and
  * non-standard element-wise vector operation can be performed using a reduction/transformation
  * operator.  As long as the individual sub-vectors are large enough, reduction/transformation
  * operators will be nearly as efficient as specialized operations for most vector subclasses.
  * Similarly to \c VectorWithOp::apply_reduction(), the \c apply_transformation() method allows
  * clients to include only subsets of elements in a reduction/transformation operation.
  *
  * In addition to being able to create non-mutable (\c const) abstract sub-views of a vector
  * object thorugh the \c VectorWithOp interface, this interface allows the creation of
  * mutable (non-<tt>const</tt>) sub-views using \c sub_view().  Also, in addition to being
  * able to extract an explicit non-mutable view of some (small?) sub-set of elements, this
  * interface allows a client to either extract a explicit mutable sub-views using
  * \c get_sub_vector() or to set sub-vectors using \c set_sub_vector(). As much
  * as possible, abstract views should be preferred (i.e. \c sub_view()) over explict views (i.e.
  * get_sub_vector() and set_sub_vector()).
  *
  * There are only three pure virtual methods that a concreate \c VectorWithOpMutable
  * subclass must override.  The \c space() and \c apply_reduction() methods from the \c VectorWithOp
  * base class inteface must be defined.  Also, as mentioned above, the \c apply_transforamtion()
  * method must be overridden and defined.
  *
  * The non-mutable (<tt>const</tt>) <tt>sub_view(...)</tt> method from the <tt>VectorWithOp</tt>
  * interface has a default implementation defined here that will be adequate for most subclasses.
  */
class VectorWithOpMutable
	: virtual public VectorBaseMutable
	, virtual public VectorWithOp
{
public:

	///
	using VectorWithOp::get_sub_vector;
	///
	using VectorWithOp::free_sub_vector;

	///
	typedef MemMngPack::ref_count_ptr<VectorWithOpMutable>    vec_mut_ptr_t;

	/** @name Pure virtual methods (must be overridden by subclass) */
	//@{

	///
	/** Apply a reduction/transformation,operation over a set of vectors:
	 * <tt>op(op(v[0]...v[nv-1],(*this),z[0]...z[nz-1]),(*reduct_obj)) -> (*this),z[0]...z[nz-1],(*reduct_obj)</tt>.
	 *
	 * The first mutable vector in the argument list to <tt>op</tt> will be
	 * <tt>this</tt> vector then followed by those in <tt>targ_vecs[k]</tt>, <tt>k = 0...num_targ_vecs-1</tt>
	 * Therefore, the number of mutable vectors passed to <tt>op\ref RTOpPack::RTOp::apply_op ".apply_op(...)"</tt>
	 * will be <tt>num_targ_vecs+1</tt>.
	 *
	 * See <tt>\ref VectorWithOp::apply_reduction "apply_reduction(...)" for a discussion of the significance of the
	 * arguments <tt>first_ele</tt>, <tt>sub_dim</tt> and <tt>global_offset</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> See <tt>\ref VectorWithOp::apply_reduction "apply_reduction(...)".
	 * </ul>
	 *
	 * @param  op	[in] Reduction/transformation operator to apply over each sub-vector
	 *				and assemble the intermediate targets into <tt>reduct_obj</tt> (if
	 *              <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
	 * @param  num_vecs
	 *				[in] Number of nonmutable vectors in <tt>vecs[]</tt>.  If <tt>vecs==NULL</tt>
	 *				then this argument is ignored but should be set to zero.
	 * @param  vecs
	 *				[in] Array (length <tt>num_vecs</tt>) of a set of pointers to
	 *				nonmutable vectors to include in the operation.
	 *				The order of these vectors is significant to <tt>op</tt>.  if <tt>vecs==NULL</tt>,
	 *              then <tt>op.apply_op(...)</tt> is called with no non-mutable sub-vector arguments.
	 * @param  num_targ_vecs
	 *				[in] Number of mutable vectors in <tt>targ_vecs[]</tt>.  If <tt>targ_vecs==NULL</tt>
	 *				then this argument is ignored but should be set to zero.
	 * @param  targ_vecs
	 *				[in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
	 *				mutable vectors to include in the operation. The order of these vectors
	 * 		        is significant to <tt>op</tt>.  If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with
	 *				only one mutable vector (<tt>*this</tt>).
	 * @param  reduct_obj
	 *				[in/out] Target object of the reduction operation.
	 *				This object must have been created by the <tt>op.reduct_obj_create_raw(&reduct_obj)</tt>
	 *              function first.  The reduction operation will be added to <tt>(*reduct_obj)</tt> if
	 *              <tt>(*reduct_obj)</tt> has already been through a reduction.  By allowing the info in
	 *              <tt>(*reduct_obj)</tt> to be added to the reduction over all of these vectors, the reduction
	 *              operation can be accumulated over a set of abstract vectors	which can be useful for implementing
	 *              composite vectors for instance.  If <tt>op.get_reduct_type_num_entries(...)</tt> returns
	 *              <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
	 *              <tt>reduct_obj</tt> should be set to #RTOp_REDUCT_OBJ_NULL and no reduction will be performed.
	 * @param  first_ele
	 *				[in] (default = 1) The index of the first element in <tt>this</tt> to be included.
	 * @param  sub_dim
	 *              [in] (default = 0) The number of elements in these vectors to include in the reduction/transformation
	 *              operation.  The value of <tt>sub_dim == 0</tt> means to include all available elements.
	 * @param  global_offset
	 *				[in] (default = 0) The offset applied to the included vector elements.
	 */
	virtual void apply_transformation(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele = 1, const index_type sub_dim = 0, const index_type global_offset = 0
		) = 0;
	
	//@}

	/** @name Virtual methods with default implementations */
	//@{

	///
	/** Assign the elements of <tt>this</tt> vector to a scalar.
	 *
	 * The default implementation of this function uses a transforamtion operator class
	 * (see RTOp_TOp_assign_scalar.h) and calls this->apply_transformation() .
	 */
	virtual VectorWithOpMutable& operator=(value_type alpha);

	///
	/** Assign the elements of a vector to <tt>this</tt>.
	 *
	 * The default implementation of this function uses a transforamtion operator class
	 * (see RTOp_TOp_assign_vectors.h) and calls <tt>this</tt>->apply_transformation() .
	 */
	virtual VectorWithOpMutable& operator=(const VectorWithOp& v);

	///
	/** Default implementation calls <tt>operator=((const &VectorWithOp)v)</tt>.
	 */
	virtual VectorWithOpMutable& operator=(const VectorWithOpMutable& v);

	///
	/** Set a specific element of a vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>1 <= i <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get(i) == val</tt>
	 * </ul>
	 *
	 * The default implementation uses a transforamtion operator
	 * class (see RTOp_TOp_set_ele.h) and calls <tt>this</tt>->apply_transforamtion().
	 *
	 * @param  i    [in] Index of the element value to set.
	 * @param  val  [in] Value of the element to set.
	 */
	virtual void set_ele( index_type i, value_type val );

	///
	/** Create a mutable abstract view of a vector object.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released by <tt>ref_count_ptr<></tt>.  This function is declared as
	 * non-constant because the object returned has the capacity to alter <tt>this</tt>
	 * object.
	 *
	 * The compatibility of sub-views goes along with the compatibility of sub-spaces
	 * (see VectorSpace).  For example, given the vector objects where
	 * <tt>x.space().is_compatible(y.space()) == true</tt> then if
	 * <tt>x.space().sub_space(rng1)->is_compatible(*y.space().sub_space(rng2)) == true</tt>
	 * then the sub-vector views <tt>*x.sub_view(rng1)</tt> and <tt>*y.sub_view(rng2)</tt>
	 * should be compatible and can be combined in vector operations.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>rng.in_range(this->dim()) == true</tt> (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.  It is allowed for the vector implementation to refuse to
	 * create arbitrary views in which case this function will return
	 * <tt>return.get() == NULL</tt>. In most applications, only specific views are
	 * every required.  The default implementation uses the subclass VectorWithOpSubView
	 * to represent any arbitrary sub-view but this can be inefficient if the sub-view is very
	 * small compared this this full vector space but not necessarily.  Note that the underlying
	 * vector <tt>this</tt> is not guarrenteed to show the changes made the sub-view
	 * <tt>*return .get()</tt> until the smart reference counted pointer <tt>return</tt> is
	 * destroyed.
	 */
	virtual vec_mut_ptr_t sub_view( const Range1D& rng );

	///
	/** Inline member function that simply calls <tt>this->sub_view(Range1D(l,u))</tt>.
	 */
	vec_mut_ptr_t sub_view( const index_type& l, const index_type& u );

	///
	/** Create a clone of this vector objet.
	 *
	 * The vector object returned in a smart reference counted pointer to a functional copy of
	 * the current vector object.  The vector object <tt>this</tt> and the vector returned by
	 * this method can be modified independently.
	 *
	 * The default implementation of this function calls on <tt>this->space().create_member()</tt> and
	 * then copies over the elements from <tt>this</tt> using <tt>operator=()</tt>.
	 */
	virtual vec_mut_ptr_t clone() const;

	///
	/** Get a mutable explicit view of a sub-vector.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released with a call to \c release_sub_vector().
	 *
	 * Note that calling this operation might require some internal
	 * allocations and temporary memory.  Therefore, it is critical
	 * that <tt>this->release_sub_vector(sub_vec)</tt> is called to
	 * clean up memory and avoid memory leaks after the sub-vector
	 * is used.
	 *
	 * If <tt>this->get_sub_vector(...,sub_vec)</tt> was previously
	 * called on <tt>sub_vec</tt> then it may be possible to reuse this
	 * memory if it is sufficiently sized.  The user is
	 * encouraged to make multiple calls to <tt>this->get_sub_vector(...,sub_vec)</tt>
	 * before <tt>this->release_sub_vector(sub_vec)</tt> to finally
	 * clean up all of the memory.  Of course the same <tt>sub_vec</tt> object must be
	 * passed to the same vector object for this to work correctly.
	 *
	 * Changes to the underlying sub-vector are not guarrenteed to become permanent
	 * until <tt>this->get_sub_vector(...,sub_vec)</tt> or <tt>this->commit_sub_vector(sub_vec)</tt>
	 * is called.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
	 *      (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * This method has a default implementation based on a vector reduction operator
	 * class (see RTOp_ROp_get_sub_vector.h) and calls ::apply_reduction<tt>(...)</tt>.
	 * Note that the footprint of the reduction object (both internal and external state)
	 * will be O(<tt>rng.size()</tt>).  For serial applications this is faily adequate and will
	 * not be a major performance penalty.  For parallel applications, this will be
	 * a terrible implementation and must be overridden if <tt>rng.size()</tt> is large at all.
	 * If a subclass does override this method, it must also override <tt>release_sub_vector()</tt>
	 * which has a default implementation which is a companion to this method's default
	 * implementation.
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sub_vec  [in/out] Mutable view of the sub-vector.  Prior to the
	 *                  first call <tt>RTOp_mutable_sub_vector_null(sub_vec)</tt> must
	 *                  have been called for the correct behavior.  Technically
	 *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
	 *                  only by calling <tt>this->commit_sub_vector(sub_vec)</tt>.
	 */
	virtual void get_sub_vector( const Range1D& rng, RTOp_MutableSubVector* sub_vec );

	///
	/** Free a mutable explicit view of a sub-vector.
	 *
	 * The sub-vector view must have been allocated by \c this->get_sub_vector() first.
	 *
	 * This method has a default implementation which is a companion to the default implementation
	 * for <tt>get_sub_vector(...)</tt>.  If <tt>get_sub_vector(...)</tt> is overridden by a subclass then
	 * this method must be overridden also!
	 *
	 *	@param	sub_vec
	 *				[in/out] The memory refered to by <tt>sub_vec->values</tt>
	 *				and <tt>sub_vec->indices</tt> will be released if it was allocated
	 *				and <tt>*sub_vec</tt> will be zeroed out using
	 *				<tt>RTOp_mutable_sub_vector_null(sub_vec)</tt>.
	 */
	virtual void commit_sub_vector( RTOp_MutableSubVector* sub_vec );

	///
	/** Set a specific sub-vector.
	 *
	 * After this function returns, the corresponding elements in <tt>this</tt> vector object will be
	 * set equal to those in the input vector (the post conditions are obvious).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>sub_vec.global_offset + sub_dim <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * The default implementation of this operation uses a transformation operator class
	 * (see RTOp_TOp_set_sub_vector.h) and calls <tt>apply_transforamtion()</tt>.  Be forewarned
	 * however, that the operator objects state data (both internal and external) will be 
	 * O(<tt>sub_vec.sub_nz</tt>).  For serial applications, this is entirely adequate.  For parallel
	 * applications this will be very bad!
	 *
	 * @param  sub_vec  [in] Represents the elements in the subvector to be set.
	 */
	virtual void set_sub_vector( const RTOp_SubVector& sub_vec );

	//@}

	/** @name Overridden from VectorWithOp */
	//@{

	///
	/** Default implementation calls <tt>this->sub_view()</tt> (non-<tt>const</tt>) and then
	 * performs an cast to <tt>vec_ptr_t</tt>.
	 *
	 * This function override is actually needed here for another reason.  Without, the
	 * override, the non-const version defined in this interface hides the const version
	 * defined in VectorWithOp.
	 */
	vec_ptr_t sub_view( const Range1D& rng ) const;

	//@}

	/** @name Overrriden from VectorBaseMutable */
	//@{

	///
	/** Calls <tt>operator=(0.0)</tt>.
	 */
	void zero();

	///
	/** Calls <tt>this->apply_transformation()</tt> with an operator class
	 * (see RTOp_TOp_axpy.h).
	 */
	void axpy( value_type alpha, const VectorBase& x );

	//@}

}; // end class VectorWithOpMutable

// ////////////////////////////////////////////////
// Inline members

inline
VectorWithOpMutable::vec_mut_ptr_t
VectorWithOpMutable::sub_view( const index_type& l, const index_type& u )
{
	return this->sub_view(Range1D(l,u));
}

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_MUTABLE_H
