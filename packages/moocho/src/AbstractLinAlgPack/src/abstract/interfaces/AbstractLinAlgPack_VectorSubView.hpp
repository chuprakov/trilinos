// ///////////////////////////////////////////////////////////////////
// VectorWithOpSubView.h

#ifndef VECTOR_WITH_OP_SUB_VIEW_H
#define VECTOR_WITH_OP_SUB_VIEW_H

#include "VectorWithOp.h"
#include "VectorSpaceSubSpace.h"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a defualt sub-view implementation for a MatrixWithOp
 * object.
 *
 * Not all of the methods from VectorWithOp are overridden, only those that
 * need to be or may result in better performance.
 *
 * The default constructor and copy constructors are allowed but the default
 * assignment operator is not allowed.
 *
 * ToDo: Finish documentation!
 */
class VectorWithOpSubView : virtual public VectorWithOp {
public:

	///
	/** Calls #this->initialize(...)#.
	 */
	VectorWithOpSubView( const vec_ptr_t& vec, const Range1D& rng );

	///
	/** Initialize.
	 *
	 * Constructs a view of the vector this = vec(rng).
	 *
	 * @param  vec    [in] The original full vector.
	 * @param  rng    [in] The range of element that #this# vector will represent.
	 */
	void initialize( const vec_ptr_t& vec, const Range1D& rng );

	///
	const VectorSpaceSubSpace& space_impl() const;

	/** @name Overridden from VectorWithOp */
	//@{

	///
	const VectorSpace& space() const;
	///
	index_type dim() const;
	///
	void apply_reduction(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type global_offset, const index_type sub_dim
		) const;
	///
	value_type get_ele(index_type i) const;
	///
	vec_ptr_t sub_view( const Range1D& rng ) const;
	///
	void get_sub_vector(
		const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const;
	///
	void free_sub_vector( RTOp_SubVector* sub_vec ) const;

	//@}

private:

	vec_ptr_t                  full_vec_;   // If full_vec_.get() == NULL, the vector is uninitalized (dim == 0)
	VectorSpaceSubSpace        space_;      // The space that this vector belongs to

	// Not defined and not to be called
	VectorWithOpSubView& operator=(const VectorWithOpSubView&);
	
}; // end class VectorWithOpSubView

// /////////////////////////////////////////////
// Inline members

inline
const VectorSpaceSubSpace& VectorWithOpSubView::space_impl() const
{
	return space_;
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_SUB_VIEW_H
