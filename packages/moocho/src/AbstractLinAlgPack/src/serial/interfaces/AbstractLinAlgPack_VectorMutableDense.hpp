// /////////////////////////////////////////////////////////////////
// VectorWithOpMutableDense.h

#ifndef VECTOR_WITH_OP_MUTABLE_DENSE_H
#define VECTOR_WITH_OP_MUTABLE_DENSE_H

#include "VectorSpaceDense.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "LinAlgPack/include/VectorClass.h"
#include "ref_count_ptr.h"
#include "ReleaseResource.h"

namespace SparseLinAlgPack {

///
/** Vector "Adaptor" subclass for <tt>LinAlgPack::VectorSlice</tt> objects.
 *
 * ToDo: Finish Documentation!
 */
class VectorWithOpMutableDense : public VectorWithOpMutable {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		ResourceManagementPack::ReleaseResource>  release_resource_ptr_t;

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorWithOpMutableDense(
		VectorSlice                        v
		,const release_resource_ptr_t&     v_release
		);

	///
	/** Initialize with a Dense vector.
	 */
	void initialize(
		VectorSlice                        v
		,const release_resource_ptr_t&     v_release
		);

	///
	VectorSlice vec();
	///
	const VectorSlice vec() const;
	///
	release_resource_ptr_t& vec_release();
	///
	const release_resource_ptr_t& vec_release() const;

	/** @name Overriddenn from VectorWithOp */
	//@{

	///
	const VectorSpace& space() const;
	///
	void apply_reduction(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const VectorWithOp** vecs
		,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		) const;
	///
	index_type dim() const;
	///
	value_type get_ele(index_type i) const;
	///
	void get_sub_vector(
		const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const;
	///
	void free_sub_vector( RTOp_SubVector* sub_vec ) const;

	//@}

	/** @name Overriddenn from VectorWithOpMutable */
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
	VectorWithOpMutable& operator=(value_type alpha);
	///
	VectorWithOpMutable& operator=(const VectorWithOp& v);
	///
	VectorWithOpMutable& operator=(const VectorWithOpMutable& v);
	///
	void set_ele( index_type i, value_type val );
	///
	void get_sub_vector(
		const Range1D& rng, RTOp_MutableSubVector* sub_vec );
	///
	void free_sub_vector( RTOp_MutableSubVector* sub_vec );
	///
	void set_sub_vector( const RTOp_SubVector& sub_vec );

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	VectorSlice               v_;
	release_resource_ptr_t    v_release_;
	VectorSpaceDense          space_;

	// ///////////////////////////////////////
	// Private member functions

	///
	/** Performs a reduction/transformation operation
	 */
	void apply_op(
		const RTOpPack::RTOp& op
		,const size_t num_vecs,      const VectorWithOp**       vecs
		,const size_t num_targ_vecs, VectorWithOpMutable**      targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		) const;

	// Not defined and not to be called
	VectorWithOpMutableDense(const VectorWithOpMutableDense&);
	VectorWithOpMutableDense& operator=(const VectorWithOpMutableDense&);

}; // end class VectorWithOpMutableDense

// //////////////////////////////////////
// Inline members

inline
VectorSlice
VectorWithOpMutableDense::vec()
{
	return v_;
}

inline
const VectorSlice
VectorWithOpMutableDense::vec() const
{
	return v_;
}

inline
VectorWithOpMutableDense::release_resource_ptr_t&
VectorWithOpMutableDense::vec_release()
{
	return v_release_;
}

inline
const VectorWithOpMutableDense::release_resource_ptr_t&
VectorWithOpMutableDense::vec_release() const
{
	return v_release_;
}

} // end namespace SparseLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_DENSE_H
