// ////////////////////////////////////////////////////////////
// VectorWithOpMutableSubView.h

#ifndef VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
#define VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H

#include "VectorWithOpMutable.h"
#include "VectorWithOpSubView.h"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a sub-view of a VectorWithOpMutable object.
 *
 * ToDo: Finish documentation! 
 */
class VectorWithOpMutableSubView
	: virtual public VectorWithOpMutable
	, virtual public VectorWithOpSubView
{
public:

	///
	/** Calls #this->initialize(...)#.
	 */
	VectorWithOpMutableSubView( const vec_mut_ptr_t& vec, const Range1D& rng );

	///
	/** Initialize.
	 *
	 * Constructs a view of the vector this = vec(rng).
	 *
	 * @param  vec    [in] The original full vector.
	 * @param  rng    [in] The range of element that #this# vector will represent.
	 */
	void initialize( const vec_mut_ptr_t& vec, const Range1D& rng );

	/** @name Overridden from VectorWithOp */
	//@{

	/// Overridden to pick VectorWithOpSubView::sub_view(rng).
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
		,const index_type global_offset, const index_type sub_dim
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

}; // end class VectorWithOpMutableSubView

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_SUB_VIEW_H
