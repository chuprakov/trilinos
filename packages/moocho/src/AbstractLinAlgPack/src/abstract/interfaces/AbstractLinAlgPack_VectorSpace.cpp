// //////////////////////////////////////////////////////////////////////
// VectorSpace.cpp
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

#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorSpaceSubSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/InnerProductDot.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpace::VectorSpace( const inner_prod_ptr_t& inner_prod )
{
	this->inner_prod(inner_prod);
}

void VectorSpace::inner_prod( const inner_prod_ptr_t& inner_prod )
{
	if(inner_prod.get()) {
		inner_prod_ = inner_prod;
	} else {
		inner_prod_ = MemMngPack::rcp(new InnerProductDot());
	}
}

const VectorSpace::inner_prod_ptr_t
VectorSpace::inner_prod() const
{
	return inner_prod_;
}

// Virtual functions with default implementations

VectorSpace::vec_mut_ptr_t
VectorSpace::create_member(const value_type& alpha) const
{
	namespace mmp = MemMngPack;
	vec_mut_ptr_t vec = this->create_member();
	*vec = alpha;
	return vec;
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpace::create_members(size_type num_vecs) const
{
	return MemMngPack::null;
}

VectorSpace::space_ptr_t
VectorSpace::sub_space(const Range1D& rng_in) const
{
	namespace mmp = MemMngPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > dim, std::out_of_range
		,"VectorSpace::sub_space(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return space_ptr_t( this, false );
	return mmp::rcp(
		new VectorSpaceSubSpace(
			mmp::rcp( this, false )
			,rng ) );
}

// Overridden from AbstractFactory<>

VectorSpace::obj_ptr_t VectorSpace::create() const
{
	return this->create_member();
}

} // end namespace AbstractLinAlgPack
