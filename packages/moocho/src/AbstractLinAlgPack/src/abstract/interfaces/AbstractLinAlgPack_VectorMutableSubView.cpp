// ////////////////////////////////////////////////////////////////
// VectorWithOpMutableSubView.cpp
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

#include "AbstractLinAlgPack/include/VectorWithOpMutableSubView.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

VectorWithOpMutableSubView::VectorWithOpMutableSubView( const vec_mut_ptr_t& vec, const Range1D& rng )
{
	this->initialize(vec,rng);
}

void VectorWithOpMutableSubView::initialize( const vec_mut_ptr_t& vec, const Range1D& rng )
{
	namespace rcp = MemMngPack;
	VectorWithOpSubView::initialize(vec,rng);
	full_vec_ = vec;
	this->has_changed();
}

void VectorWithOpMutableSubView::set_uninitialized()
{
	VectorWithOpSubView::set_uninitialized();
	full_vec_ = MemMngPack::null;
	this->has_changed();
}

// Overriddend form VectorWithOp

VectorWithOp::vec_ptr_t VectorWithOpMutableSubView::sub_view( const Range1D& rng ) const
{
	return VectorWithOpSubView::sub_view(rng); // Had to override to resolve conflicit!
}

// Overridden from VectorWithOpMutable

void VectorWithOpMutableSubView::apply_transformation(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const VectorWithOp** vecs
	,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	)
{
#ifdef _DEBUG
	const index_type this_dim = this->dim();
	THROW_EXCEPTION(
		sub_dim_in < 0
		|| !(1 <= first_ele_in && first_ele_in <= this_dim)
		|| ( sub_dim_in > 0 && (sub_dim_in - (first_ele_in - 1) > this_dim) )
		, std::logic_error
		,"VectorWithOpMutableSubView::apply_transformation(...): Error, first_ele_in = "
		<< first_ele_in << ", global_offset_in = " << global_offset_in
		<< ", sub_dim_in = " << sub_dim_in << " and this->dim() = this_dim  = "
		<< this_dim << " are not compatible." );
#endif
	const index_type this_offset = space_impl().rng().lbound() - 1;
	const index_type
		this_sub_dim = ( sub_dim_in 
						 ? sub_dim_in
						 : space_impl().rng().size() - (first_ele_in - 1) );
	full_vec_->apply_transformation(
		op, num_vecs, vecs, num_targ_vecs, targ_vecs, reduct_obj
		,(index_type)(this_offset + first_ele_in)     // first_ele
		,this_sub_dim                                 // sub_dim
		,global_offset_in                             // global_dim
		);
}

void VectorWithOpMutableSubView::set_ele( index_type i, value_type val )
{
	space_impl().validate_range(Range1D(i,i));
	full_vec_->set_ele( space_impl().rng().lbound() + i - 1, val );
}

VectorWithOpMutable::vec_mut_ptr_t
VectorWithOpMutableSubView::sub_view( const Range1D& rng_in )
{
	namespace rcp = MemMngPack;
	const size_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
	space_impl().validate_range(rng);
	if( rng.lbound() == 1 && rng.ubound() == this_dim )
		return rcp::rcp(this,false); // Do not own memory!
	const index_type this_offset = space_impl().rng().lbound() - 1;
	return rcp::rcp(
		new VectorWithOpMutableSubView(
			full_vec_
			,Range1D( 
				this_offset  + rng.lbound()
				,this_offset + rng.ubound() )
			) );
}

void VectorWithOpMutableSubView::set_sub_vector( const RTOp_SubVector& sub_vec_in )
{
	const index_type  this_offset = space_impl().rng().lbound() - 1;
	RTOp_SubVector    sub_vec = sub_vec_in;
	sub_vec.global_offset += this_offset;
	full_vec_->set_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
