// //////////////////////////////////////////
// VectorWithOp.h

#ifndef VECTOR_WITH_OP_H
#define VECTOR_WITH_OP_H

#include "VectorBase.h"
#include "LinAlgPack/include/Range1D.h"

namespace AbstractLinAlgPack {

class VectorWithOpMutable;

///
/** Abstract interface for immutable vectors {abstract}.
  *
  * This interface contains a mimimal set of operations.  The main feature
  * of this interface is the operation \Ref{apply_reduction}#(...)#.
  * Almost every standard (i.e. BLAS) and non-standard operation that
  * can be performed on a set of vectors without changing (mutating)
  * the vectors can be performed through reduction operators.  Standard
  * vector operations could be included in this interface and allow
  * for specialized implementations but in general, assuming the
  * sub-vectors are large enough, such implementations
  * would not be significantly faster than those implemented through
  * reduction/transformation operators.  Therefore, the only operation that an
  * immutable abstract vector must implement is the #apply_reduction(...)#
  * method.
  */
class VectorWithOp : virtual public VectorBase {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorWithOp>   vec_ptr_t;

	///
	/** Fetch an element in the vector {abstract}.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #1 <= i <= this->dim()# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  i  [in]  Index of the element value to get.
	 */
	virtual RTOp_value_type get_ele(RTOp_index_type i) const = 0;

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
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.
	 */
	virtual vec_ptr_t create_sub_view( const LinAlgPack::Range1D& rng ) const = 0;


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
	 * is used in a reduction or transformation operation.
	 *
	 * If #this->extract_sub_vector(...,sub_vec)# was previously
	 * called on #sub_vec# then it may be possible to reuse this
	 * memory if it is sufficiently sized.  The user is
	 * encouraged to make multiple calls to #this->extract_sub_vector(...,sub_vec)#
	 * before #this->release_sub_vector(sub_vec)# to finally
	 * clean up all of the memory.  Of course the same #sub_vec# object must be
	 * passed to the same vector object for this to work correctly.
	 *
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sub_vec  [in/out] View of the sub-vector.  Prior to the
	 *                  first call \Ref{RTOp_sub_vector_null}#(sub_vec)# must
	 *                  have been called for the correct behavior.
	 */
	virtual void get_sub_vector( const LinAlgPack::Range1D& rng, RTOp_SubVector* sub_vec ) const = 0;

	///
	/** Release an explicit view of a sub-vector {abstract}.
	 *
	 * The sub-vector view must have been allocated by
	 * this->\Ref{extract_sub_vector}#(...,sub_vec)# first.
	 *
	 *	@param	sub_vec
	 *				[in/out] The memory refered to by #sub_vec->values#
	 *				and #sub_vec->indices# will be released if it was allocated
	 *				and #*sub_vec# will be zeroed out using
	 *				\Ref{RTOp_sub_vector_null}#(sub_vec)#.
	 */
	virtual void release_sub_vector( RTOp_SubVector* sub_vec ) const = 0;
	
	//@}

};

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_H
