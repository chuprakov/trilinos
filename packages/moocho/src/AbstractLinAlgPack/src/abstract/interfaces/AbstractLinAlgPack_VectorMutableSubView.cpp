// ////////////////////////////////////////////////////////////////
// VectorMutableSubView.cpp
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

#include "VectorMutableSubView.hpp"
#include "Teuchos_TestForException.hpp"
#include "WorkspacePack.hpp"
#include "dynamic_cast_verbose.hpp"

namespace AbstractLinAlgPack {

VectorMutableSubView::VectorMutableSubView( const vec_mut_ptr_t& vec, const Range1D& rng )
{
	this->initialize(vec,rng);
}

void VectorMutableSubView::initialize( const vec_mut_ptr_t& vec, const Range1D& rng )
{
	namespace rcp = MemMngPack;
	VectorSubView::initialize(vec,rng);
	full_vec_ = vec;
	this->has_changed();
}

void VectorMutableSubView::set_uninitialized()
{
	VectorSubView::set_uninitialized();
	full_vec_ = Teuchos::null;
	this->has_changed();
}

// Overriddend form Vector

Vector::vec_ptr_t VectorMutableSubView::sub_view( const Range1D& rng ) const
{
	return VectorSubView::sub_view(rng); // Had to override to resolve conflicit!
}

// Overridden from VectorMutable

void VectorMutableSubView::set_ele( index_type i, value_type val )
{
	space_impl().validate_range(Range1D(i,i));
	full_vec_->set_ele( space_impl().rng().lbound() + i - 1, val );
}

VectorMutable::vec_mut_ptr_t
VectorMutableSubView::sub_view( const Range1D& rng_in )
{
	namespace rcp = MemMngPack;
	const size_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
	space_impl().validate_range(rng);
	if( rng.lbound() == 1 && rng.ubound() == this_dim )
		return Teuchos::rcp(this,false); // Do not own memory!
	const index_type this_offset = space_impl().rng().lbound() - 1;
	return Teuchos::rcp(
		new VectorMutableSubView(
			full_vec_
			,Range1D( 
				this_offset  + rng.lbound()
				,this_offset + rng.ubound() )
			) );
}

void VectorMutableSubView::set_sub_vector( const RTOpPack::SparseSubVector& sub_vec_in )
{
	const index_type            this_offset = space_impl().rng().lbound() - 1;
	RTOpPack::SparseSubVector   sub_vec = sub_vec_in;
	sub_vec.setGlobalOffset( sub_vec.globalOffset() + this_offset );
	full_vec_->set_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
