// //////////////////////////////////////////
// VectorWithOp.hpp
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

#ifndef VECTOR_WITH_OP_H
#define VECTOR_WITH_OP_H

#include <iosfwd>

#include "AbstractLinAlgPackTypes.hpp"
#include "RTOpPack/src/RTOpCpp.hpp"
#include "Range1D.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract interface for immutable, finite dimensional, coordinate vectors {abstract}.
 *
 * This interface contains a mimimal set of operations.  The main feature
 * of this interface is the operation \c apply_reduction().
 * Almost every standard (i.e. BLAS) and non-standard element-wise operation that
 * can be performed on a set of coordinate vectors without changing (mutating)
 * the vectors can be performed through reduction operators.  More standard
 * vector operations could be included in this interface and allow
 * for specialized implementations but in general, assuming the
 * sub-vectors are large enough, such implementations
 * would not be significantly faster than those implemented through
 * reduction/transformation operators.  There are some operations however
 * that can not always be efficiently with reduction/transforamtion operators
 * and a few of these important methods are included in this interface.  The
 * <tt>apply_reduction()</tt> method allows to client to specify a sub-set
 * of the vector elements to include in reduction/transformation operation.
 * This greatly increases the generality of this vector interface as vector
 * objects can be used as sub objects in larger composite vectors and sub-views
 * of a vector can be created.
 *
 * This interface allows clients to create sub-views of a vector using \c sub_view()
 * that in turn are fully functional <tt>%VectorWithOp</tt> objects.  This functionality
 * is supported by default by using a default vector subclass \c VectorWithOpSubView which
 * in turn calls <tt>apply_reduction()</tt> but the client need not ever worry about
 * how this is done.
 *
 * This interface also allows a client to extract a sub-set of elements in an
 * explicit form as an \c RTOp_SubVector object using the method \c get_sub_vector().
 * In general, this is very bad thing to do and should be avoided at all costs.
 * However, there are some applications where this is needed and therefore it is
 * supported.  The default implementation of this method uses a reduction/transformation
 * operator with <tt>apply_reduction()</tt> in order to extract the needed elements.
 *
 * In order to create a concreate subclass of this interface, only two
 * methods must be overridden.  The \c space() method must be overridden which in turn
 * requires defining a concreate \c VectorSpace class (which has only two pure virtual
 * methods).  And, as mentioned above, the \c apply_reduction() method must be overridden
 * as well.
 *
 * The fact that this interface defines \c space() which returns a \c VectorSpace object
 * (which in turn can create mutable vectors) implies that for every possible vector object,
 * it is possible to associate with it a mutable vector object that can be the target
 * of transformation operations.  This is not a serious limitation.  For any
 * application area, mutable vectors should be able to defined and should be
 * usable with the non-mutable vectors.
 *
 * This interface includes methods for the common vector norms: \c norm_1(),
 * \c norm_2(), \c norm_inf().  The default implementation of this class uses reduction
 * operator classes (See RTOp_ROp_norms.h) and caches away the values of the
 * norms that are computed since it is common that the norms will be accessed many
 * times before a vector is changed.  The operations in any subclass that modifies
 * the underlying vector must call the method <tt>this-></tt>has_changed() in order
 * to alert this implementation that the norms are no longer valid.
 *
 * ToDo: Add example code!
 */
class VectorWithOp {
public:

	///
	typedef MemMngPack::ref_count_ptr<const VectorWithOp>   vec_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<VectorWithOpMutable>  vec_mut_ptr_t;
	///
	enum ESparseOrDense { SPARSE, DENSE };

	///
	VectorWithOp();
	///
	virtual ~VectorWithOp() {}

	/** @name Pure virtual methods (must be overridden by subclass) */
	//@{

	///
	/** Return the vector space that this vector belongs to.
	 *
	 * Note that the vectors space object returned is specifically bound to this
	 * vector object.  The vector space object returned should only be considered
	 * to be transient and may become invalid if <tt>this</tt> is modified in some way
	 * (though some other interface).
	 */
	virtual const VectorSpace& space() const = 0;

	///
	/** Apply a reduction/transformation,operation over a set of vectors:
	 * <tt>op(op((*this),v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
	 *
	 * The first nonmutable vector in the argument list to <tt>op</tt> will be <tt>this</tt>
	 * vector then followed by those in <tt>vecs[k], k = 0...num_vecs-1</tt>.  Therefore, the number
	 * of nonmutable vectors passed to <tt>op\ref RTOpPack::RTOp::apply_op ".apply_op(...)"</tt> will be
	 * <tt>num_vecs+1</tt>.
	 *
	 * The vector to be represented <tt>v</tt> by <tt>this</tt> and passed to the
	 * <tt>op\ref RTOpPack::RTOp::apply_op ".apply_op(...)"</tt> method is:
	 \verbatim

	 v(k + global_offset) = this->get_ele(first_ele + k - 1)
         , for k = 1 ... sub_dim
     \endverbatim
     * The other vector arguments are represented similarly.  The situation where
	 * <tt>first_ele == 1</tt> and <tt>global_offset > 1</tt> corresponds to the
	 * case where the vectors are representing consitituent vectors in a larger
	 * aggregrate vector.  The situation where <tt>first_ele > 1</tt> and
	 * <tt>global_offset == 0</tt> is for when a sub-view of the vectors are being
	 * treated as full vectors.  Other combinations of these arguments is
	 * possible.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>num_vecs > 0</tt>] <tt>vecs[k]->space().is_compatible(this->space()) == true</tt>
	 *          , for <tt>k = 0...num_vecs-1</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> [<tt>num_targ_vecs > 0</tt>] <tt>targ_vecs[k]->space().is_compatible(this->space()) == true</tt>
	 *          , for <tt>k = 0...num_targ_vecs-1</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>1 <= first_ele <= this->dim()</tt> (throw <tt>std::out_of_range</tt>)
	 * <li> <tt>global_offset >= 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>sub_dim - (first_ele - 1) <= this->dim()</tt> (throw <tt>std::length_error</tt>).
	 * </ul>
	 *
	 * @param  op	[in] Reduction/transformation operator to apply over each sub-vector
	 *				and assemble the intermediate targets into <tt>reduct_obj</tt> (if
	 *              <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
	 * @param  num_vecs
	 *				[in] Number of nonmutable vectors in <tt>vecs[]</tt>.
	 *              If <tt>vecs==NULL</tt> then this argument is ignored but should be set to zero.
	 * @param  vecs
	 *				[in] Array (length <tt>num_vecs</tt>) of a set of pointers to
	 *				nonmutable vectors to include in the operation.
	 *				The order of these vectors is significant to <tt>op</tt>.
	 *				If <tt>vecs==NULL</tt> then <tt>op</tt> is called with the
	 *				single vector represented by <tt>this</tt> object.
	 * @param  num_targ_vecs
	 *				[in] Number of mutable vectors in <tt>targ_vecs[]</tt>.
	 *              If <tt>targ_vecs==NULL</tt>	then this argument is ignored but should be set to zero.
	 * @param  targ_vecs
	 *				[in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
	 *				mutable vectors to include in the operation.
	 *				The order of these vectors is significant to <tt>op</tt>.
	 *				If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with no mutable vectors.
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
	 *
	 * <b> Note the subclass implementors </b>
	 *
	 * It is imporatant that all transformed vectors have the method \c has_changed() called on them
	 * durring the implementation of this function.  This can be done by calling
	 * \c finalize_apply_reduction(...) just before the \c apply_reduction() method returns.
	 * For example, the implementation for a vector subclass <tt>MyVector</tt> would look like:
	 \code
	 
	 void MyVector::apply_reduction(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		) const
	 {
	     ...
		 finalize_apply_reduction(num_targ_vecs,targ_vecs);
	 }
	 \endcode
	 */
	virtual void apply_reduction(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele = 1, const index_type sub_dim = 0, const index_type global_offset = 0
		) const = 0;

	//@}

	/** @name Miscellaneous virtual methods with default implementations */
	//@{

	///
	/** Return the dimension of this vector.
	 *
	 * It is allowed for a vector to return a dimension of <tt>0</tt> in which case
	 * the vector should be considered uninitialized in which the client should
	 * not call any of the member functions (except space()).  The default implementation
	 * returns <tt>this->space().dim()</tt>.
	 */
	virtual index_type dim() const;

	///
	/** Return the number of nonzero elements in the vector.
	 *
	 * The default implementation just uses a reduction operator
	 * with the apply_reduction<tt>(...)</tt> method (See
	 * RTOp_ROp_num_nonzeros.h).
	 */
	virtual index_type nz() const;

	///
	/** Virtual output function.
	  *
	  * The default implementation just uses get_sub_vector<tt>(...)</tt> to convert to
	  * a dense vector and then prints this.
	  *
	  * ToDo: Finish documentation!
	  *
	  * @param  out        [in/out] Receives the output.
	  * @param  print_dim  [in] (default = true) If true, then the dimension is printed
	  *                    first followed by a newline.
	  * @param  newline    [in] (default = true) If true, then a newline is printed after
	  *                    the last element is printed.
	  * @param  global_offset
	  *                    [in] (default = 0) The offset added to the vector element
	  *                    indexes when they are printed.
	  */
	virtual std::ostream& output(
		std::ostream& out, bool print_dim = true, bool newline = true
		,index_type global_offset = 0 ) const;

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
	/** Fetch an element in the vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>1 <= i <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * The default implementation uses a C reduction operator class
	 * (See RTOp_ROp_get_ele.h C).
	 *
	 * @param  i  [in]  Index of the element value to get.
	 */
	virtual value_type get_ele(index_type i) const;

	///
	/** Create an abstract view of a vector object .
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released by <tt>ref_count_ptr<></tt>.
	 *
	 * It is important to understand what the minimum postconditions are for the sub vector objects
	 * returned from this method.  If two vector objects <tt>x</tt> and <tt>y</tt> are compatible (possibly of
	 * different types) it is assumed that <tt>*x.sub_view(rng)</tt> and <tt>*y.sub_view(rng)</tt>
	 * will also be compatible vector objects no mater what range <tt>rng</tt> represents.  However,
	 * if <tt>i1 < i2 < i3 < i4</tt> with <tt>i2-i1 == i4-i3</tt>, then in general, one can not expect
	 * that the vector objects <tt>*x.sub_view(Range1D(i2,i1))</tt> and
	 * <tt>*x.sub_view(Range1D(i4,i5))</tt> will be compatible objects.  For some vector
	 * implementaions they may be (i.e. serial vectors) but for others they most
	 * certainly will not be (i.e. parallel vectors).  This limitation must be kept in
	 * mind by all vector subclass implementors and vector interface clients.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
	 *      (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>[return.get() != NULL] return->get_ele(i) == this->get_ele(i+rng.lbound()-1)</tt>
	 *       , for <tt>i = 1...rng.size()</tt>.
	 * </ul>
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.  It
	 *              is allowed for <tt>rng.full_range() == true</tt> in which case it implicitly
	 *              treated as <tt>rng = [1,this->dim()]</tt>.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.  It is allowed for subclasses to return  <tt>return->get() == NULL</tt>
	 * for some selections
	 * of <tt>rng</tt>.  Only some <tt>rng</tt> ranges may be allowed but they will be appropriate for the
	 * application at hand.  However, a very good implementation should be able to
	 * accommodate any valid <tt>rng</tt> that meets the basic preconditions.  The default
	 * implementation uses the subclass \c VectorWithOpSubView to represent any arbitrary
	 * sub-view but this can be inefficient if the sub-view is very small compared this this
	 * full vector space but not necessarily.
	 */
	virtual vec_ptr_t sub_view( const Range1D& rng ) const;

	//@}

	///
	/** Inline member function that simply calls <tt>this->sub_view(Range1D(l,u))</tt>.
	 */
	vec_ptr_t sub_view( const index_type& l, const index_type& u ) const;

	/** @name Vector norms */
	//@{

	///
	/** One norm. <tt>||v||_1 = sum( |v(i)|, i = 1,,,this->dim() )</tt>
	 */
	virtual value_type norm_1() const;
	///
	/** Two norm. <tt>||v||_2 = sqrt( sum( v(i)^2, i = 1,,,this->dim() ) )</tt>
	 */
	virtual value_type norm_2() const;
	///
	/** Infinity norm.  <tt>||v||_inf = max( |v(i)|, i = 1,,,this->dim() )</tt>
	 */
	virtual value_type norm_inf() const;
	
	//@}

	/** @name Inner product */
	//@{

	///
	/** Return the inner product of <tt>*this</tt> with <tt>v</tt>.
	 *
	 * @return Returns <tt>this->space().inner_prod()->inner_prod(*this,v)</tt>
	 */
	virtual value_type inner_product( const VectorWithOp& v ) const;

	//@}

	/** @name Explicit sub-vector access */
	//@{

	///
	/** Get a non-mutable explicit view of a sub-vector.
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
	 * If a subclass does override this method, it must also override <tt>release_sub_vector(...)</tt>
	 * which has a default implementation which is a companion to this method's default
	 * implementation.
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sparse_or_dense
	 *                  [in] If <tt>SPARSE</tt> then <tt>*sub_vec</tt> will be a sparse vector 
	 *                  (i.e. <tt>sub_vec->indices != NULL</tt>) if there are any zero elements.
	 *                  If <tt>DENSE</tt> then the retured <tt>*sub_vec</tt> will be dense (i.e.
	 *                  <tt>sub_vec->indices == NULL</tt>) and will contain all of the zero
	 *                  as well as the nonzero elements.
	 * @param  sub_vec  [in/out] View of the sub-vector.  Prior to the
	 *                  first call <tt>RTOp_sub_vector_null(sub_vec)</tt> must
	 *                  have been called for the correct behavior.  Technically
	 *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
	 *                  only by calling <tt>this->free_sub_vector(sub_vec)</tt>.
	 */
	virtual void get_sub_vector(
		const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const;

	///
	/** Free an explicit view of a sub-vector.
	 *
	 * The sub-vector view must have been allocated by this->get_sub_vector() first.
	 *
	 * This method has a default implementation which is a companion to the default implementation
	 * for <tt>get_sub_vector(...)</tt>.  If <tt>get_sub_vector(...)</tt> is overridden by a subclass then
	 * this method must be overridden also!
	 *
	 *	@param	sub_vec
	 *				[in/out] The memory refered to by <tt>sub_vec->values</tt>
	 *				and <tt>sub_vec->indices</tt> will be released if it was allocated
	 *				and <tt>*sub_vec</tt> will be zeroed out using
	 *				<tt>RTOp_sub_vector_null(sub_vec)</tt>.
	 */
	virtual void free_sub_vector( RTOp_SubVector* sub_vec ) const;

	//@}

	///
	/** Must be called by any vector subclass that modifies this vector
	 * object!
	 *
	 * The way to use this method by subclasses is to call it when ever
	 * there is a chance that the vector may have changed.  Therefore, this
	 * method should be called in every non-const member function in every
	 * subclass.  This is a little bit of a pain but overall this is a good
	 * design in that it allows for efficient cacheing of information for
	 * multiple retreval.  For example, if the subclass <tt>SomeVector</tt> has cashed
	 * data and has a method <tt>SomeVector::foo()</tt> may modify the
	 * vector then <tt>SomeVector</tt> should override the method <tt>has_changed()</tt> and its
	 * implementation should look someting likde like this!
	 \verbatim
	 void SomeVector::has_changed()
	 {
	     BaseClass::has_changed(); // Called on most direct subclass that
		                           // has overridden this method as well.
	    ...  // Reinitialize your own cached data to uninitialized!
	 }
	 \endverbatim
	 */
	virtual void has_changed() const;

protected:

	/** @name Protected helper functions */
	//@{

	///
	/** This method usually needs to be called by subclasses at the
	 * end of the \c apply_reduction() method implementation to
	 * insure that \c has_changed() is called on the transformed
	 * vector objects.
	 */
	virtual void finalize_apply_reduction(
		const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		) const;

	//@}

private:

	mutable index_type  num_nonzeros_;  ///< < 0 == not initialized, > 0 == already calculated
	mutable value_type  norm_1_, norm_2_, norm_inf_;  ///< < 0 == not initialized, > 0 == already calculated

}; // end class MatrixWithOp

// ////////////////////////////////////////////////
// Inline members

inline
VectorWithOp::vec_ptr_t
VectorWithOp::sub_view( const index_type& l, const index_type& u ) const
{
	return this->sub_view(Range1D(l,u));
}

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_H
