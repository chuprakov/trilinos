// /////////////////////////////////////////////////////////////////////////
// VectorWithOpSubView.cpp
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

#include <assert.h>

#include <stdexcept>

#include "AbstractLinAlgPack/include/VectorWithOpMutableSubView.h"
#include "WorkspacePack.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace AbstractLinAlgPack {

VectorWithOpSubView::VectorWithOpSubView( const vec_ptr_t& vec, const Range1D& rng )
{
	initialize(vec,rng);
}

void VectorWithOpSubView::initialize( const vec_ptr_t& vec, const Range1D& rng )
{
	namespace rcp = MemMngPack;
#ifdef _DEBUG
	THROW_EXCEPTION(
		vec.get() == NULL, std::invalid_argument
		,"VectorWithOpSubView::initialize(...) : Error!" );
#endif
	typedef VectorSpace::space_ptr_t   space_ptr_t;
	space_.initialize( rcp::rcp(&vec->space(),false), rng );
	full_vec_ = vec;
	this->has_changed();
}

void VectorWithOpSubView::set_uninitialized()
{
	space_.set_uninitialized();
	full_vec_ = MemMngPack::null;
	this->has_changed();
}

// Overridden from VectorWithOp

const VectorSpace& VectorWithOpSubView::space() const
{
	return space_;
}

index_type VectorWithOpSubView::dim() const
{
	return space_.dim();
}

void VectorWithOpSubView::apply_reduction(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const VectorWithOp** vecs
	,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	int k;
	const index_type this_dim = this->dim();
#ifdef _DEBUG
	THROW_EXCEPTION(
		sub_dim_in < 0
		|| !(1 <= first_ele_in && first_ele_in <= this_dim)
		|| ( sub_dim_in > 0 && (sub_dim_in - (first_ele_in - 1) > this_dim) )
		, std::logic_error
		,"VectorWithOpSubView::apply_reduction(...): Error, first_ele_in = "
		<< first_ele_in << ", global_offset_in = " << global_offset_in
		<< ", sub_dim_in = " << sub_dim_in << " and this->dim() = this_dim  = "
		<< this_dim << " are not compatible." );
#endif
	const index_type this_offset = space_impl().rng().lbound() - 1;
	const index_type
		this_sub_dim = ( sub_dim_in 
						 ? sub_dim_in
						 : space_impl().rng().size() - (first_ele_in - 1)
			           );
	wsp::Workspace<const VectorWithOp*>    vecs_full(wss,num_vecs);
	for( k = 0; k < num_vecs; ++k )
		vecs_full[k] = dyn_cast<const VectorWithOpSubView>(*vecs[k]).full_vec().get();
	wsp::Workspace<VectorWithOpMutable*>   targ_vecs_full(wss,num_targ_vecs);
	for( k = 0; k < num_targ_vecs; ++k )
		targ_vecs_full[k] = dyn_cast<VectorWithOpMutableSubView>(*targ_vecs[k]).full_vec().get();
	full_vec_->apply_reduction(
		op
		,num_vecs,      num_vecs      ? &vecs_full[0]      : NULL
		,num_targ_vecs, num_targ_vecs ? &targ_vecs_full[0] : NULL
		,reduct_obj
		,this_offset + first_ele_in     // first_ele
		,this_sub_dim                   // sub_dim
		,global_offset_in               // global_dim
		);
}

value_type VectorWithOpSubView::get_ele(index_type i) const
{
	space_.validate_range(Range1D(i,i));
	return full_vec_->get_ele( space_.rng().lbound() + i - 1 );
}

VectorWithOp::vec_ptr_t
VectorWithOpSubView::sub_view( const Range1D& rng_in ) const
{
	namespace rcp = MemMngPack;
	const index_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range(rng_in,1,this_dim);
	space_.validate_range(rng);
	const index_type this_offset = space_.rng().lbound() - 1;
	return rcp::rcp(
		new VectorWithOpSubView(
			full_vec_
			,Range1D( 
				this_offset  + rng.lbound()
				,this_offset + rng.ubound() )
			) );
}

void VectorWithOpSubView::get_sub_vector(
	const Range1D& rng_in, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		!sub_vec, std::logic_error
		,"VectorWithOpSubView::get_sub_vector(...): Error!" ) ;
#endif
	const index_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range(rng_in,1,this_dim);
	space_.validate_range(rng);
	const index_type this_offset = space_.rng().lbound() - 1;
	full_vec_->get_sub_vector( rng + this_offset, sparse_or_dense, sub_vec );
	sub_vec->global_offset -= this_offset;
}

void VectorWithOpSubView::free_sub_vector( RTOp_SubVector* sub_vec ) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		!sub_vec, std::logic_error
		,"VectorWithOpSubView::free_sub_vector(...): Error!" ) ;
#endif
	const index_type this_offset = space_.rng().lbound() - 1;
	sub_vec->global_offset += this_offset;
	full_vec_->free_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
