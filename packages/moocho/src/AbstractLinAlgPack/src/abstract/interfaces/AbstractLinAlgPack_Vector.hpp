// //////////////////////////////////////////
// VectorWithOp.h

#ifndef VECTOR_WITH_OP_H
#define VECTOR_WITH_OP_H

#include "VectorBase.h"

namespace AbstractLinAlgPack {

class VectorWithOpMutable;

///
/** Abstract interface for immutable coordinate vectors {abstract}.
  *
  * This interface contains a mimimal set of operations.  The main feature
  * of this interface is the operation \Ref{apply_reduction}#(...)#.
  * Almost every standard (i.e. BLAS) and non-standard operation that
  * can be performed on a set of coordinate vectors without changing (mutating)
  * the vectors can be performed through reduction operators.  More standard
  * vector operations could be included in this interface and allow
  * for specialized implementations but in general, assuming the
  * sub-vectors are large enough, such implementations
  * would not be significantly faster than those implemented through
  * reduction/transformation operators.  There are some operations however
  * that can not always be efficiently with reduction/transforamtion operators
  * and a few of these important methods are added to this interface.
  *
  * In order to create a concreate subclass of this interface, only three
  * methods absolutely must be overridden.  The method \Ref{dim}#()#
  * gives the dimension of the vector at hand. As mentioned above,
  * the \Ref{apply_reduction}#(...)# method must be
  * overridden.  The last method that must be overridden is \Ref{create_sub_view}#(...)#.
  * This method allows clients to perform operations on sub-vectors within the
  * whole vector.  This makes the #apply_reduction(...)# method even more powerful
  * but makes vector implementations slightly more complicated.
  *
  * The default implementation of this class caches away the values of the
  * norms that are computed.  The operations in any subclass that modifies
  * the underlying vector must call the method #this->has_changed()# in order
  * to alert this implementation that the norms are no longer valid.
  */
class VectorWithOp : virtual public VectorBase {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorWithOp>   vec_ptr_t;

	///
	enum ESparseOrDense { SPARSE, DENSE };

	///
	VectorWithOp();

	/** @name Pure virtual methods (must be overridden by subclass).
	 */
	//@{

	/// Return the dimension of this vector
	virtual RTOp_index_type dim() const = 0;

	///
	/** Apply a reduction/transformation,operation over a set of vectors:
	 * #op(op((*this),v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)#.
	 *
	 * The first nonmutable vector in the argument list to #op# will be
	 * #this# vector then followed by those in #vecs[k]#, k = 0...#num_vecs-1
	 * Therefore, the number of nonmutable vectors passed to
	 * \Ref{RTOp_apply_op}#(...)# will be #num_vecs+1#.
	 *
	 * Preconditions:\begin{itemize}
	 * \item [#num_vecs > 0#] #vecs[k]->dim() == this->dim()#, for k = 0...num_vecs-1
	 * \item [#num_targ_vecs > 0#] #vecs[k]->dim() == this->dim()#, for k = 0...num_targ_vecs-1
	 * \end{itemize}
	 *
	 * @param  op	[in] Reduction operator to apply over each sub-vector
	 *				and assemble the intermediate targets into #reduct_obj#.
	 *	@param  num_vecs
	 *				[in] Number of nonmutable vectors in #vecs[]#.  If #vecs==NULL#
	 *				then this argument is ignored but should be set to zero.
	 *	@param  vecs
	 *				[in] Array (length #num_vecs#) of a set of pointers to
	 *				nonmutable vectors to include in the operation.
	 *				The order of these vectors is significant to #op#.
	 *				If #vecs==NULL# then #op# is called with the
	 *				single vector represented by #this# object.
	 *	@param  num_targ_vecs
	 *				[in] Number of mutable vectors in #targ_vecs[]#.  If #targ_vecs==NULL#
	 *				then this argument is ignored but should be set to zero.
	 *	@param  targ_vecs
	 *				[in] Array (length #num_targ_vecs#) of a set of pointers to
	 *				mutable vectors to include in the operation.
	 *				The order of these vectors is significant to #op#.
	 *				If #targ_vecs==NULL# then #op# is called with no mutable vectors.
	 *	@param  reduct_obj
	 *				[in/out] Target object of the reduction operation.
	 *				This object must have been created by the
	 *				#op#.\Ref{reduct_obj_create_raw}#(&reduct_obj)# function
	 *				first.  The reduction operation will be added
	 *				to #(*reduct_obj)# if #(*reduct_obj)# has already been through a
	 *				reduction.  By allowing the info in #(*reduct_obj)#
	 *				to be added to the reduction over all of these
	 *				vectors, the reduction operation can
	 *				be accumulated over a set of abstract vectors
	 *				which can be useful for implementing concatenated
	 *				vectors for instance.
	 *             If #op#.\Ref{get_reduct_type_num_entries}#(...)# returns
	 *             #num_values == 0#, #num_indexes == 0# and #num_chars == 0#
	 *             then #reduct_obj# should be set to #RTOp_REDUCT_OBJ_NULL#
	 *             and no reduction will be performed.	  
	 */
	virtual void apply_reduction( const RTOpPack::RTOp& op
		, const size_t num_vecs, const VectorWithOp** vecs
		, const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		, RTOp_ReductTarget reduct_obj ) const = 0;

	///
	/** Create an abstract view of a vector object {abstract}.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released by #ref_count_ptr<>#.
	 *
	 * It is important to understand what the minimum postconditions are for the sub vector objects
	 * returned from this method.  If two vector objects #x# and #y# are compatible (possibly of
	 * different types) it is assumed that #*x.create_sub_view(rng)# and #*y.create_sub_view(rng)#
	 * will also be compatible vector objects no mater what range #rng# represents.  However,
	 * if #i1 < i2 < i3 < i4# with #i2-i1 == i4-i3#, then in general, one can not expect
	 * that the vector objects #*x.create_sub_view(Range1D(i2,i1))# and
	 * #*x.create_sub_view(Range1D(i4,i5))# will be compatible objects.  For some vector
	 * implementaions they may be (i.e. serial vectors) but for others they most
	 * certainly will not be (i.e. parallel vectors).  This limitation must be kept in
	 * mind by all vector subclass implementors and vector interface clients.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.  It is allowed that #return->get() == NULL# for some selections
	 * of #rng#.  Only some #rng# ranges may be allowed but they will be appropriate for the
	 * application at hand.  However, a very good implementation should be able to
	 * accommodate any valid #rng# that meets the basic preconditions.
	 */
	virtual vec_ptr_t create_sub_view( const Range1D& rng ) const = 0;

	//@}

	/** @name Virtual methods with default implementations based on
	 * reduction/transforamtion operators and \Ref{apply_reduction}#(...)#.
	 */
	//@{

	///
	/** Fetch an element in the vector {abstract}.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #1 <= i <= this->dim()# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * The default implementation uses a C reduction operator class
	 * (See the RTOp_ROp_get_ele C operator class).
	 *
	 * @param  i  [in]  Index of the element value to get.
	 */
	virtual RTOp_value_type get_ele(RTOp_index_type i) const;

	/** @name Vector norms
	 *
	 * These member functions have default implementations based on
	 * reduction operator classes.  The default implementation of this
	 * class caches the value of the norms since it is common that the
	 * norms will be accessed many times before a vector is changed.
	 */

	///
	/** One norm.
	 * @memo #||v||_1 = sum( |v(i)|, i = 1,,,this->dim() )#
	 */
	virtual RTOp_value_type norm_1() const;
	///
	/** Two norm.
	 * @memo #||v||_2 = sqrt( sum( v(i)^2, i = 1,,,this->dim() ) )#
	 */
	virtual RTOp_value_type norm_2() const;
	///
	/** Infinity norm.
	 * @memo #||v||_inf = max( |v(i)|, i = 1,,,this->dim() )#
	 */
	virtual RTOp_value_type norm_inf() const;
	
	//@}

	/** @name Explicit sub-vector access.
	 *
	 * These member functions can be used to extract a explicit view
	 * of any subvector in the overall vector.  Note that this may be
	 * a very bad thing to do with many vector subclasses (i.e. parallel
	 * and out-of-core vectors).  Allowing a user to create an explict view
	 * of the elements allows great flexibility but must be used with care
	 * and only when absolutly needed.
	 */
	//@{

	///
	/** Get an explicit view of a sub-vector {abstract}.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released with a call to \Ref{release_sub_vector}#(...)#.
	 *
	 * Note that calling this operation might require some internal
	 * allocations and temporary memory.  Therefore, it is critical
	 * that this->\Ref{release_sub_vector}#(sub_vec)# is called to
	 * clean up memory and avoid memory leaks after the sub-vector
	 * is used.
	 *
	 * If #this->get_sub_vector(...,sub_vec)# was previously
	 * called on #sub_vec# then it may be possible to reuse this
	 * memory if it is sufficiently sized.  The user is
	 * encouraged to make multiple calls to #this->get_sub_vector(...,sub_vec)#
	 * before #this->release_sub_vector(sub_vec)# to finally
	 * clean up all of the memory.  Of course the same #sub_vec# object must be
	 * passed to the same vector object for this to work correctly.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * This method has a default implementation based on a vector reduction operator
	 * class (see RTOp_ROp_get_sub_vector) and calls \Ref{apply_reduction}#(...)#.
	 * Note that the footprint of the reduction object (both internal and external state)
	 * will be O(#rng.size()#).  For serial applications this is faily adequate and will
	 * not be a major performance penalty.  For parallel applications, this will be
	 * a terrible implementation and must be overridden if #rng.size()# is large at all.
	 * If a subclass does override this method, it must also override \Ref{release_sub_vector}#(...)#
	 * which has a default implementation which is a companion to this method's default
	 * implementation.
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sparse_or_dense
	 *                  [in] If #SPARSE# then #*sub_vec# will be a sparse vector 
	 *                  (i.e. #sub_vec->indices != NULL#) if there are any zero elements.
	 *                  If #DENSE# then the retured #*sub_vec# will be dense (i.e.
	 *                  #sub_vec->indices == NULL#) and will contain all of the zero
	 *                  as well as the nonzero elements.
	 * @param  sub_vec  [in/out] View of the sub-vector.  Prior to the
	 *                  first call \Ref{RTOp_sub_vector_null}#(sub_vec)# must
	 *                  have been called for the correct behavior.  Technically
	 *                  #*sub_vec# owns the memory but this memory can be freed
	 *                  only by calling #this->free_sub_vector(sub_vec)#.
	 */
	virtual void get_sub_vector(
		const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const;

	///
	/** Free an explicit view of a sub-vector {abstract}.
	 *
	 * The sub-vector view must have been allocated by
	 * this->\Ref{get_sub_vector}#(...,sub_vec)# first.
	 *
	 * This method has a default implementation which is a companion to the default implementation
	 * for \Ref{get_sub_vector}#(...)#.  If #get_sub_vector(...)# is overridden by a subclass then
	 * this method must be overridden also!
	 *
	 *	@param	sub_vec
	 *				[in/out] The memory refered to by #sub_vec->values#
	 *				and #sub_vec->indices# will be released if it was allocated
	 *				and #*sub_vec# will be zeroed out using
	 *				\Ref{RTOp_sub_vector_null}#(sub_vec)#.
	 */
	virtual void free_sub_vector( RTOp_SubVector* sub_vec ) const;

	///
	/** Must be called by any vector subclass that modifies this vector
	 * object!
	 *
	 * The way to use this method by subclass is to call it when ever
	 * there is a chance that the vector may have changed.  Therefore, this
	 * method should be called in every non-const member function in every
	 * subclass.  This is a little bit of a pain but overall this is a good
	 * design.  For example, if the subclass #SomeVector# has cashed
	 * data as well and some method #SomeVector::foo()# may modify the
	 * vector then #SomeVector# should override this method and its
	 * implementation should look like this!
	 \begin{verbatim}
	 void SomeVector::has_changed()
	 {
	     BaseClass::has_changed(); // Called on most direct subclass that
		                           // has overridden this method as well.
	    ...  // Reinitialize your own cached data to uninitialized!
	 }
	 \end{verbatim}
	 */
	virtual void has_changed() const;

	//@}

	/** @name Overridden from \Ref{VectorBase} */
	//@{
	
	///
	/** Calls #apply_reduction(...)# with an operator class object
	 * of type #RTOp_ROp_dot_prod#.
	 */
	RTOp_value_type inner_product(  const VectorBase& ) const;

	//@}

private:

	mutable RTOp_value_type  norm_1_, norm_2_, norm_inf_;   // < 0, not initialized, > 0 already calculated

};

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_H
