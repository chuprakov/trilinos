// //////////////////////////////////////////
// VectorWithOpMutable.h

#ifndef VECTOR_WITH_OP_MUTABLE_H
#define VECTOR_WITH_OP_MUTABLE_H

#include "VectorWithOp.h"

namespace AbstractLinAlgPack {

///
/** Abstract interface for objects that represent a space for vectors.
 *
 */
class VectorSpace {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorSpace>    space_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<VectorWithOpMutable>  vec_mut_ptr_t;

	///
	virtual ~VectorSpace() {}

	///
	/** Return the dimmension of the vector space.
	 */
	virtual RTOp_index_type dim() const = 0;

	///
	/** Compare the compatibility of two vector spaces.
	 *
	 * If this function returns true, then vectors created from
	 * either of the vector spaces will be compatible and can
	 * be combined in vector operations.
	 */
	virtual bool is_compatible(const VectorSpace& another_vector_space) const = 0;

	///
	/** Create a member of the vector space.
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector object.  On construction #return->get_ele(i)# will
	 * be equal to #0.0# for #i = 1...this->dim()#.  Note that #&return->space()#
	 * does not have to be equal to #this#.  
	 */
	virtual vec_mut_ptr_t create_member() const = 0;

	///
	/** Create a subspace of the current vector space.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng  [in] The range of the elements to extract a vector sub-space.
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector space object.  Note that the vector object returned
	 * by #this->sub_space(rng).create_member()# should be exactly equivalent
	 * to the vector returned by
	 * #this->create_member()->create_sub_view(rng)->space()->create_member()#.
	 * It is allowed for the implementation to return #return->get() == NULL#
	 * for arbitrary values of #rng#.  Only some #rng# ranges may be allowed
	 * but they will be appropriate for the application at hand.
	 */
	virtual space_ptr_t sub_space(const LinAlgPack::Range1D& rng) const = 0;

}; // end class VectorSpace

///
/** Abstract interface for mutable vectors {abstract}.
  *
  * Objects of this type can act as a target vector of
  * reduction/transformation operations.  Similarly to 
  * \Ref{VectorWithOp} this interface contains a single
  * method \Ref{apply_transformation}#(...)# that allows
  * users to apply user defined transformation operators.
  * Every standard (i.e. BLAS) and non-standard element-wise
  * vector operation can be performed using a transformation
  * operator and therefore the operation #apply_transformation(...)#
  * is the only method an mutable abstract vector must implement.
  * As long as the individual sub-vectors are large enough,
  * transformation operators will be nearly as efficient as
  * specialized operations for most vector subclasses so why
  * even include any other methods?
  */
class VectorWithOpMutable : virtual public VectorWithOp {
public:

	///
	typedef VectorSpace::vec_mut_ptr_t    vec_mut_ptr_t;

	///
	/** Set a specific element of a vector.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #1 <= i <= this->dim()# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  i    [in] Index of the element value to set.
	 * @param  val  [in] Value of the element to set.
	 */
	virtual void set_ele( RTOp_index_type i, RTOp_value_type val ) = 0;

	///
	/** Apply a reduction/transformation,operation over a set of vectors:
	  * #op(op(v[0]...v[nv-1],(*this),z[0]...z[nz-1]),(*reduct_obj)) -> (*this),z[0]...z[nz-1],(*reduct_obj)#.
	  *
	  * The first mutable vector in the argument list to #op# will be
	  * #this# vector then followed by those in #targ_vecs[k]#, k = 0...#num_targ_vecs-1
	  * Therefore, the number of mutable vectors passed to
	  * \Ref{RTOp_apply_op}#(...)# will be #num_targ_vecs+1#.
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
	  *				\Ref{RTOp_reduct_obj_create}#(op,reduct_obj)# function
	  *				first.  The reduction operation will be added
	  *				to #(*reduct_obj)# if #(*reduct_obj)# has already been through a
	  *				reduction.  By allowing the info in #(*reduct_obj)#
	  *				to be added to the reduction over all of these
	  *				vectors, the reduction operation can
	  *				be accumulated over a set of abstract vectors
	  *				which can be useful for implementing concatenated
	  *				vectors for instance.
 	  *             If \Ref{RTOp_get_reduct_type_num_entries}#(op,...)# returns
	  *             #num_values == 0#, #num_indexes == 0# and #num_chars == 0#
	  *             then #reduct_obj# should be set to #RTOp_REDUCT_OBJ_NULL#
	  *             and no reduction will be performed.	  
	  */
	virtual void apply_transformation( const RTOpPack::RTOp& op
		, const size_t num_vecs, const VectorWithOp** vecs
		, const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		, RTOp_ReductTarget reduct_obj ) = 0;

	///
	/** Return the vector space that this vector belongs to {abstract}.
	 *
	 * Note that the vectors space object returned is specifically bound to this
	 * vector object.  The vector space object returned should only be considered
	 * to be transient and may become invalid if #this# is modified in some way
	 * (but not through the VectorWithOpMutable interface).
	 */
	virtual const VectorSpace& space() const = 0;
	
	///
	/** Create a clone of this vector objet {abstract}.
	 *
	 * The vector object returned in the smart reference counted pointer
	 * is a functional copy of the current vector object.  The vector object
	 * #this# and the vector returned by this method can be modified independently.
	 */
	virtual vec_mut_ptr_t clone() const = 0;

	///
	/** Create a mutable abstract view of a vector object {abstract}.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released by #ref_count_ptr<>#.  This function is declared as
	 * non-constant because the object returned has the capacity to alter #this#
	 * object.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.  It is allowed for the vector implementation to refuse to
	 * create arbitrary views in which case this function will return #return.get() == NULL#.
	 * In most applications, only specific views are every required.
	 */
	virtual vec_mut_ptr_t create_sub_view( const LinAlgPack::Range1D& rng ) = 0;

	///
	/** Set a specific sub-vector {abstract}.
	 *
	 * After this function returns, the corresponding elements in #this# vector object will be
	 * set equal to those in the input vector (the post conditions are obvious).
	 *
	 * Preconditions:\begin{itemize}
	 * \item #sub_vec.global_offset + sub_dim <= this->dim()# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  sub_vec  [in] Represents the elements in the subvector to be set.
	 */
	virtual void set_sub_vector( const RTOp_SubVector& sub_vec ) = 0;

};

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_MUTABLE_H
