// /////////////////////////////////////////////////////////////////
// VectorWithOpMutableDense.cpp
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

#include <typeinfo>
#include <stdexcept>

#include "SparseLinAlgPack/src/VectorWithOpMutableDense.hpp"
#include "SparseLinAlgPack/src/VectorDenseEncap.hpp"
#include "SparseLinAlgPack/src/GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack/src/apply_op_helper.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

#ifdef _DEBUG
#define CLASS_MEMBER_PTRS \
const VectorWithOpMutableDense   *_this = this; \
const DVectorSlice                *_v; \
const release_resource_ptr_t     *_v_release; \
const VectorSpaceSerial          *_space;
#else
#define CLASS_MEMBER_PTRS
#endif

namespace SparseLinAlgPack {

VectorWithOpMutableDense::VectorWithOpMutableDense(
	const size_type                   dim
	)
	:space_(dim)
{
	CLASS_MEMBER_PTRS
	this->initialize(dim);
}

VectorWithOpMutableDense::VectorWithOpMutableDense(
	DVectorSlice                        v
	,const release_resource_ptr_t&     v_release
	)
	:space_(v.dim())
{
	CLASS_MEMBER_PTRS
	this->initialize(v,v_release);
}

void VectorWithOpMutableDense::initialize(
	const size_type                   dim
	)
{
	CLASS_MEMBER_PTRS
	namespace rcp = MemMngPack;
	namespace rmp = MemMngPack;
	typedef rcp::ref_count_ptr<DVector> vec_ptr_t;
	vec_ptr_t vec_ptr = rcp::rcp(new DVector(dim));
	this->initialize(
		(*vec_ptr)()
		,rcp::rcp(
			new rmp::ReleaseResource_ref_count_ptr<DVector>(
				vec_ptr
				)
			)
		);
}

void VectorWithOpMutableDense::initialize(
	DVectorSlice                        v
	,const release_resource_ptr_t&     v_release
	)
{
	CLASS_MEMBER_PTRS
	v_.bind(v);
	v_release_ = v_release;
	space_.initialize(v.dim());
	this->has_changed();
}

// Overridden from VectorWithOp

const VectorSpace& VectorWithOpMutableDense::space() const
{
	CLASS_MEMBER_PTRS
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
	CLASS_MEMBER_PTRS
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
	CLASS_MEMBER_PTRS
	const size_type  this_dim = v_.dim();
	const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorWithOpMutableDense::get_sub_vector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
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
	CLASS_MEMBER_PTRS
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
	CLASS_MEMBER_PTRS
	v_ = alpha;
	this->has_changed();
	return *this;
}

VectorWithOpMutable&
VectorWithOpMutableDense::operator=(const VectorWithOp& v)
{
	CLASS_MEMBER_PTRS
	if( const VectorWithOpMutableDense *vp = dynamic_cast<const VectorWithOpMutableDense*>(&v) )
		v_ = vp->v_;
	else
		return VectorWithOpMutable::operator=(v); // Try the default implementation?
	this->has_changed();
	return *this;
}

VectorWithOpMutable&
VectorWithOpMutableDense::operator=(const VectorWithOpMutable& v)
{
	CLASS_MEMBER_PTRS
	if( const VectorWithOpMutableDense *vp = dynamic_cast<const VectorWithOpMutableDense*>(&v) )
		v_ = vp->v_;
	else
		return VectorWithOpMutable::operator=(v); // Try the default implementation?
	this->has_changed();
	return *this;
}

void VectorWithOpMutableDense::set_ele( index_type i, value_type val )
{
	CLASS_MEMBER_PTRS
	v_(i) = val;
	this->has_changed();
}

VectorWithOpMutableDense::vec_mut_ptr_t
VectorWithOpMutableDense::sub_view( const Range1D& rng_in )
{
	CLASS_MEMBER_PTRS
	namespace rcp = MemMngPack;
	const size_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorWithOpMutableDense::sub_view(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
#endif
	if( rng == Range1D(1,this_dim) )
		return rcp::rcp( this, false );
	this->has_changed(); // This will result in a change in the vector
	return rcp::rcp( new VectorWithOpMutableDense( v_(rng), rcp::null ) ); 
}

void VectorWithOpMutableDense::get_sub_vector(
	const Range1D& rng_in, RTOp_MutableSubVector* sub_vec )
{
	CLASS_MEMBER_PTRS
	const size_type  this_dim = v_.dim();
	const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorWithOpMutableDense::get_sub_vector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
#endif
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

void VectorWithOpMutableDense::commit_sub_vector( RTOp_MutableSubVector* sub_vec )
{
	CLASS_MEMBER_PTRS
	RTOp_mutable_sub_vector_null( sub_vec ); // No memory to deallocate!
	this->has_changed(); // Be aware of any final changes!
}

void VectorWithOpMutableDense::set_sub_vector( const RTOp_SubVector& sub_vec )
{
	CLASS_MEMBER_PTRS
	VectorWithOpMutable::set_sub_vector(sub_vec); // ToDo: Provide specialized implementation?
}

void VectorWithOpMutableDense::Vp_StMtV(
	value_type                       alpha
	,const GenPermMatrixSlice        &P
	,BLAS_Cpp::Transp                P_trans
	,const VectorWithOp              &x
	,value_type                      beta
	)
{
	CLASS_MEMBER_PTRS
	VectorDenseEncap  x_de(x);
	SparseLinAlgPack::Vp_StMtV( &v_, alpha, P, P_trans, x_de(), beta );
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
	CLASS_MEMBER_PTRS
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
