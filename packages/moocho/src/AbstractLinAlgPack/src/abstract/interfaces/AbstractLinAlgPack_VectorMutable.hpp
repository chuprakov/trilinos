// //////////////////////////////////////////
// VectorWithOpMutable.h

#ifndef VECTOR_WITH_OP_MUTABLE_H
#define VECTOR_WITH_OP_MUTABLE_H

#include "VectorBaseMutable.h"
#include "VectorWithOp.h"

namespace AbstractLinAlgPack {

///
/** Abstract interface for objects that represent a space for mutable coordinate vectors.
 *
 * This class is primary an "Abstract Factory" interface.  A #VectorSpace# object can exist
 * independent from any individual #VectorWithOpMutable# object.  Or, a #VectorSpace# object
 * can have a lifetime that is dependent on a single #VectorWithOpMutable# object.  The
 * same interface can serve both roles.  
 */
class VectorSpace : public virtual VectorSpaceBase {
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
	/** Create a member of the vector space.
	 *
	 * Postconditions:\begin{itemize}
	 * \item #return->dim() == this->dim()#
	 * \item #return->get_ele(i) == 0.0#, for #i = 1...this->dim()#
	 * \item #return->space().is_compatible(*this) == true#
	 * \end{itemize}
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
	 * Postconditions:\begin{itemize}
	 * \item #return->dim() == rng->size()#
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
	 * but they will be appropriate for the application at hand.  However, a
	 * very good implementation should be able to accomidate any valid #rng#
	 * that meets the basic preconditions.
	 *
	 * Note that if two vector space objects #X# and #Y# are compatible (i.e.
	 * #X.is_compatible(Y) == true#, then it is also expected that
	 * #X.sub_space(rng)->is_compatible(*Y.sub_space(rng))# will also be true.
	 * However, in general it can not be expected that
	 * #X.sub_space(rng1)->is_compatible(*X.sub_space(rng2))#, with
	 * #rng1.size() == rng2.size()#, will be true.  For serial vectors, it may
	 * be but for parallel vectors it will most certainly not be.  Therefore, in
	 * general, don't assume that arbitrary subsets of the vector spaces will be
	 * comatible, even if the sizes of these subspaces are the same.
	 */
	virtual space_ptr_t sub_space(const Range1D& rng) const = 0;
	
	/** @name Overrriden from \Ref{VectorSpaceBase} */
	//@{
	
	///
	/** This implementaion just calls #create_member()# and then converts
	 * the smart pointer.
	 */
	vec_ptr_t new_member() const;

	//@}

}; // end class VectorSpace

///
/** Abstract interface for mutable coordinate vectors {abstract}.
  *
  * Objects of this type can act as a target vector of
  * a transformation operation.  Similarly to \Ref{VectorWithOp}
  * this interface contains very few pure virtual methods that must
  * be overridden.  However, more efficient implementations will
  * choose to override more methods.
  * 
  * The most important method is \Ref{apply_transformation}#(...)#
  * that allows clients to apply user defined reduction/transformation
  * operators.  Every standard (i.e. BLAS) and non-standard element-wise
  * vector operation can be performed using a reduction/transformation
  * operator.  As long as the individual sub-vectors are large enough,
  * transformation operators will be nearly as efficient as
  * specialized operations for most vector subclasses so why
  * even include any other methods?
  *
  * There are only few pure virtual methods that a concreate subclass
  * must override.  The #dim()# method from the \Ref{VectorWithOp} base
  * class interface must be overridden to give the dimmension of the vector.
  * The \Ref{apply_reduction}#(...)# method from the \Ref{VectorWithOp}
  * base class inteface must be defined.  The mutable (non-const)
  * \Ref{create_sub_view}#(...)# method from this interface must be
  * overriden.  The #space()# method must also be overridden which in
  * turn requires defining a concreate #VectorSpace# class.  Therefore,
  * defining a concreate #VectorWithOpMutable# subclass is significanly
  * harder that definign a concreate #VectorWithOp# class.  However,
  * note only four vector methods and three #VectorSpace# methods must
  * be overridden to define a very powerful coordinate vector class.
  *
  * The non-mutable (const) \Ref{create_sub_view}#(...)# from the
  * #VectorWithOp# interface has a default implementation defined
  * here that will be adequate for most subclasses.
  */
class VectorWithOpMutable
	: virtual public VectorBaseMutable
	, virtual public VectorWithOp
{
public:

	///
	typedef VectorSpace::vec_mut_ptr_t    vec_mut_ptr_t;

	/** @name Pure virtual methods (must be overridden by subclass).
	 */
	//@{

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
	  *				mutable vectors to include in the operation. The order of these vectors
	  * 			is significant to #op#.  If #targ_vecs==NULL# then #op# is called with
	  *				only one mutable vector (#*this#).
	  *	@param  reduct_obj
	  *				[in/out] Target object of the reduction operation.
	  *				This object must have been created by the
	  *				#op#.\Ref{reduct_obj_create_raw}#(&reduct_obj)# function
	  *				first.  The reduction operation will be added
	  *				to #(*reduct_obj)# if #(*reduct_obj)# has already been
	  *				through a reduction.  By allowing the info in #(*reduct_obj)#
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
	/** Create a mutable abstract view of a vector object {abstract}.
	 *
	 * This is only a transient view of a sub-vector that is to be immediately used
	 * and then released by #ref_count_ptr<>#.  This function is declared as
	 * non-constant because the object returned has the capacity to alter #this#
	 * object.
	 *
	 * The compatibility of sub-views goes along with the compatibility of subspaces
	 * (see \Ref{VectorSpace}).  For example, given the vector objects where
	 * #x.space().is_compatible(y.space()) == true# then if
	 * #x.space().sub_space(rng1)->is_compatible(*y.space().sub_space(rng2)) == true#
	 * then the sub-vector views #*x.create_sub_view(rng1)# and #*y.create_sub_view(rng2)#
	 * should be compatible and can be combined in vector operations.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #rng.in_range(this->dim()) == true# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * @param  rng  [in] The range of the elements to extract the sub-vector view.
	 * 
	 * @return  Returns a smart reference counted pointer to a view of the requested
	 * vector elements.  It is allowed for the vector implementation to refuse to
	 * create arbitrary views in which case this function will return
	 * #return.get() == NULL#. In most applications, only specific views are every required.
	 */
	virtual vec_mut_ptr_t create_sub_view( const Range1D& rng ) = 0;

	//@}

	/** @name Virtual methods with default implementations based on
	 * reduction/transforamtion operators and \Ref{apply_transforamtion}#(...)#.
	 */
	//@{

	///
	/** Assign the elements of this vector to a scalar.
	 *
	 * The default implementation of this function uses a transforamtion
	 * operator class (see RTOp_TOp_assign_scalar) and calls \Ref{apply_transformation}#(...)#.
	 */
	virtual VectorWithOpMutable& operator=(RTOp_value_type);

	///
	/** Assign the elements of of a vector to this.
	 *
	 * The default implementation of this function uses a transforamtion
	 * operator class (see RTOp_TOp_assign_vector) and calls \Ref{apply_transformation}#(...)#.
	 */
	virtual VectorWithOpMutable& operator=(const VectorWithOpMutable&);

	///
	/** Set a specific element of a vector.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #1 <= i <= this->dim()# (#throw std::out_of_range#)
	 * \end{itemize}
	 *
	 * The default implementation uses a transforamtion operator
	 * class (see RTOp_TOp_set_ele) and calls \Ref{apply_transforamtion}#(...)#.
	 *
	 * @param  i    [in] Index of the element value to set.
	 * @param  val  [in] Value of the element to set.
	 */
	virtual void set_ele( RTOp_index_type i, RTOp_value_type val );

	///
	/** Create a clone of this vector objet {abstract}.
	 *
	 * The vector object returned in the smart reference counted pointer
	 * is a functional copy of the current vector object.  The vector object
	 * #this# and the vector returned by this method can be modified independently.
	 *
	 * The default implementation of this function calls on #this->space().create_member()# and
	 * then copies over the elements from #this# using #operator=(...)#.
	 */
	virtual vec_mut_ptr_t clone() const;

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
	 * The default implementation of this operation uses a transformation operator class
	 * (see RTOp_TOp_set_sub_vector) and calls \Ref{apply_transforamtion}#(...)#.  Be forewarned
	 * however, that the operator objects state data (both internal and external) will be
	 * O(#sub_vec.sub_nz#).  For serial applications, this is entirely adequate.  For parallel
	 * applications this will be very bad!
	 *
	 * @param  sub_vec  [in] Represents the elements in the subvector to be set.
	 */
	virtual void set_sub_vector( const RTOp_SubVector& sub_vec );

	//@}

	/** @name Overridden from \Ref{VectorWithOp} */
	//@{

	///
	/** Default implementation calls create_sub_view (non-const) and then converts.
	 *
	 * This function override is actually needed here for another reason.  Without, the
	 * override, the non-const version defined in this interface hides the const version
	 * defined in \Ref{VectorWithOp}.
	 */
	vec_ptr_t create_sub_view( const Range1D& rng ) const;

	//@}

	/** @name Overrriden from \Ref{VectorBaseMutable} */
	//@{
	
	///
	/** Calls #this-space()#.
	 */
	const VectorSpaceBase& get_space() const;

	///
	/** Calls #operator=(0.0)#.
	 */
	void zero();

	///
	/** Calls #apply_transformation(...)# with the operator class
	 * RTOp_TOp_axpy.
	 */
	void axpy( RTOp_value_type alpha, const VectorBase& x );

	//@}

}; // end class VectorWithOpMutable

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_MUTABLE_H
