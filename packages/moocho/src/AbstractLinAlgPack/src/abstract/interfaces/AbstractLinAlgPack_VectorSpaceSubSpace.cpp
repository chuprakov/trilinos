// //////////////////////////////////////////////////////////////////////
// VectorSpaceSubSpace.cpp
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

#include "AbstractLinAlgPack_VectorSpaceSubSpace.hpp"
#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

VectorSpaceSubSpace::VectorSpaceSubSpace( const space_ptr_t& full_space, const Range1D& rng )
{
	this->initialize(full_space,rng);
}

void VectorSpaceSubSpace::initialize( const space_ptr_t& full_space, const Range1D& rng )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		full_space.get() == NULL, std::invalid_argument
		,"VectorSpaceSubSpace::initialize(...): Error!" );
#endif
	const index_type n = full_space->dim();
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!rng.full_range() && rng.ubound() > n, std::out_of_range
		,"VectorSpaceSubSpace::initialize(...): Error, "
		"rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
		"[1,vec->dim()] = [1," << n << "]" );
#endif
	full_space_ = full_space;
	rng_ = rng.full_range() ? Range1D(1,n) : rng;
}

void VectorSpaceSubSpace::set_uninitialized()
{
	full_space_ = Teuchos::null;
	rng_        = Range1D::Invalid;
}

#ifdef _DEBUG
void VectorSpaceSubSpace::validate_range(const Range1D& rng) const
{
	const index_type n = this->dim();
	TEST_FOR_EXCEPTION(
		full_space_.get() == NULL, std::logic_error
		,"VectorSpaceSubSpace::validate_range(rng): Error, Uninitialized" );
	TEST_FOR_EXCEPTION(
		full_space_.get() && !rng.full_range() && rng.ubound() > n, std::logic_error
		,"VectorSpaceSubSpace::validate_range(rng): Error, "
		"rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
		"[1,this->dim] = [1," << n << "]" );
}
#endif

// Overridden from VectorSpace

bool VectorSpaceSubSpace::is_compatible(const VectorSpace& another_space) const
{
	if( this->dim() == another_space.dim() && this->is_in_core() && another_space.is_in_core() )
		return true;
	const VectorSpaceSubSpace
		*a_space = dynamic_cast<const VectorSpaceSubSpace*>(&another_space);
	if(!a_space)
		return false;
	return
		( this->full_space_.get() == NULL && a_space->full_space_.get() == NULL )
		||
		( this->rng_ == a_space->rng_ && this->full_space_->is_compatible(*a_space->full_space_) );
}

bool VectorSpaceSubSpace::is_in_core() const
{
	return full_space_->is_in_core();
}

index_type VectorSpaceSubSpace::dim() const
{
	return full_space_.get() ? rng_.size() : 0;
}

VectorSpace::vec_mut_ptr_t VectorSpaceSubSpace::create_member() const
{
	namespace rcp = MemMngPack;
	if( full_space_.get() )
		return Teuchos::rcp(
			new VectorMutableSubView(
				full_space_->create_member(), rng_ 
				) );
	return Teuchos::null;
}

VectorSpace::space_ptr_t VectorSpaceSubSpace::clone() const
{
	namespace rcp = MemMngPack;
	if( full_space_.get() )
		return Teuchos::rcp(new VectorSpaceSubSpace( full_space_->clone(), rng_ ));
	return Teuchos::rcp(new VectorSpaceSubSpace());
}

VectorSpace::space_ptr_t VectorSpaceSubSpace::sub_space(const Range1D& rng_in) const
{
	namespace rcp = MemMngPack;
	if( full_space_.get() == NULL && rng_in == Range1D::Invalid )
		return Teuchos::rcp(this,false);
	validate_range(rng_in);
	const index_type dim         = this->dim();
	const Range1D    rng         = rng_in.full_range() ? Range1D(1,dim) : rng_in;
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return space_ptr_t( this, false );
	const index_type this_offset = rng_.lbound() - 1;
	return Teuchos::rcp(
		new VectorSpaceSubSpace(
			full_space_
			,Range1D( 
				this_offset  + rng.lbound()
				,this_offset + rng.ubound() )
			) );
}

} // end namespace AbstractLinAlgPack
