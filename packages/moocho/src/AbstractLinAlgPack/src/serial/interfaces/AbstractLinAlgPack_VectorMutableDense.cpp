// /////////////////////////////////////////////////////////////////
// VectorWithOpMutableDense.cpp

#include <typeinfo>
#include <stdexcept>

#include "SparseLinAlgPack/include/VectorWithOpMutableDense.h"
#include "AbstractLinAlgPack/include/apply_op_helper.h"
#include "ReleaseResource_ref_count_ptr.h"
#include "WorkspacePack.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

VectorWithOpMutableDense::VectorWithOpMutableDense(
	const size_type                   dim
	)
	:space_(dim)
{
	this->initialize(dim);
}

VectorWithOpMutableDense::VectorWithOpMutableDense(
	VectorSlice                        v
	,const release_resource_ptr_t&     v_release
	)
	:space_(v.dim())
{
	this->initialize(v,v_release);
}

void VectorWithOpMutableDense::initialize(
	const size_type                   dim
	)
{
	namespace rcp = ReferenceCountingPack;
	namespace rmp = ResourceManagementPack;
	typedef rcp::ref_count_ptr<Vector> vec_ptr_t;
	vec_ptr_t vec_ptr = rcp::rcp(new Vector(dim));
	this->initialize(
		(*vec_ptr)()
		,rcp::rcp(
			new rmp::ReleaseResource_ref_count_ptr<Vector>(
				vec_ptr
				)
			)
		);
}

void VectorWithOpMutableDense::initialize(
	VectorSlice                        v
	,const release_resource_ptr_t&     v_release
	)
{
	v_.bind(v);
	v_release_ = v_release;
	space_.initialize(v.dim());
	this->has_changed();
}

// Overridden from VectorWithOp

const VectorSpace& VectorWithOpMutableDense::space() const
{
	return space_;
}

void VectorWithOpMutableDense::apply_reduction(
	const RTOpPack::RTOp& op
	,const size_t num_vecs_in, const VectorWithOp** vecs_in
	,const size_t num_targ_vecs_in, VectorWithOpMutable** targ_vecs_in
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele, const index_type sub_dim, const index_type global_offset
	) const
{
 	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	wsp::Workspace<const VectorWithOp*>  vecs(wss,num_vecs_in+1);
	vecs[0] = this; // I am the first argument!
	std::copy( vecs_in, vecs_in + num_vecs_in, &vecs[1] );

	apply_op(
		op
		,num_vecs_in+1,     &vecs[0]
		,num_targ_vecs_in,  targ_vecs_in
		,reduct_obj, first_ele, sub_dim, global_offset
		);
}

index_type VectorWithOpMutableDense::dim() const
{
	return v_.dim();
}

value_type VectorWithOpMutableDense::get_ele(index_type i) const
{
	return v_(i);
}

void VectorWithOpMutableDense::get_sub_vector(
	const Range1D& rng_in, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const
{
	const size_type  this_dim = v_.dim();
	const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorWithOpMutableDense::get_sub_vector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim );
	// Just return the dense view regardless of spare_or_dense argument
	RTOp_SubVector _sub_vec;
	RTOp_sub_vector_dense(
		rng.lbound()-1                             // global_offset
		,rng.size()                                // sub_dim
		,v_.raw_ptr()+v_.stride()*(rng.lbound()-1) // values
		,v_.stride()
		,&_sub_vec
		);
	*sub_vec = _sub_vec;  // No memory has been allocated here!
}

void VectorWithOpMutableDense::free_sub_vector( RTOp_SubVector* sub_vec ) const
{
	RTOp_sub_vector_null( sub_vec ); // No memory to deallocate!
}

// Overriddenn from VectorWithOpMutable

void VectorWithOpMutableDense::apply_transformation(
	const RTOpPack::RTOp& op
	,const size_t num_vecs_in, const VectorWithOp** vecs_in
	,const size_t num_targ_vecs_in, VectorWithOpMutable** targ_vecs_in
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele, const index_type sub_dim, const index_type global_offset
	)
{
 	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	wsp::Workspace<VectorWithOpMutable*>  targ_vecs(wss,num_targ_vecs_in+1);
	targ_vecs[0] = this; // I am the first argument!
	std::copy( targ_vecs_in, targ_vecs_in + num_targ_vecs_in, &targ_vecs[1] );

	apply_op(
		op
		,num_vecs_in,         vecs_in
		,num_targ_vecs_in+1,  &targ_vecs[0]
		,reduct_obj, first_ele, sub_dim, global_offset
		);
}

VectorWithOpMutable&
VectorWithOpMutableDense::operator=(value_type alpha)
{
	v_ = alpha;
}

VectorWithOpMutable&
VectorWithOpMutableDense::operator=(const VectorWithOp& v)
{
	if( const VectorWithOpMutableDense *vp = dynamic_cast<const VectorWithOpMutableDense*>(&v) )
		v_ = vp->v_;
	else
		return VectorWithOpMutable::operator=(v); // Try the default implementation?
	return *this;
}

VectorWithOpMutable&
VectorWithOpMutableDense::operator=(const VectorWithOpMutable& v)
{
	if( const VectorWithOpMutableDense *vp = dynamic_cast<const VectorWithOpMutableDense*>(&v) )
		v_ = vp->v_;
	else
		return VectorWithOpMutable::operator=(v); // Try the default implementation?
	return *this;
}

void VectorWithOpMutableDense::set_ele( index_type i, value_type val )
{
	v_(i) = val;
}

void VectorWithOpMutableDense::get_sub_vector(
	const Range1D& rng_in, RTOp_MutableSubVector* sub_vec )
{
	const size_type  this_dim = v_.dim();
	const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorWithOpMutableDense::get_sub_vector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim );
	RTOp_MutableSubVector _sub_vec;
	RTOp_mutable_sub_vector(
		rng.lbound()-1                             // global_offset
		,rng.size()                                // sub_dim
		,v_.raw_ptr()+v_.stride()*(rng.lbound()-1) // values
		,v_.stride()
		,&_sub_vec
		);
	*sub_vec = _sub_vec;  // No memory has been allocated here!
}

void VectorWithOpMutableDense::free_sub_vector( RTOp_MutableSubVector* sub_vec )
{
	RTOp_mutable_sub_vector_null( sub_vec ); // No memory to deallocate!
}

void VectorWithOpMutableDense::set_sub_vector( const RTOp_SubVector& sub_vec )
{
	VectorWithOpMutable::set_sub_vector(sub_vec); // ToDo: Provide specialized implementation?
}

// private

void VectorWithOpMutableDense::apply_op(
	const RTOpPack::RTOp& op
	,const size_t num_vecs,      const VectorWithOp**  vecs
	,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	) const
{
	const index_type this_dim = this->dim();
#ifdef _DEBUG
	AbstractLinAlgPack::apply_op_validate_input(
		NULL, "VectorWithOpMutableDense::apply_op(...)"
		,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele_in,sub_dim_in,global_offset_in
		);
#endif
	AbstractLinAlgPack::apply_op_serial(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
		,first_ele_in,sub_dim_in,global_offset_in );
}

} // end namespace SparseLinAlgPack
