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

#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/VectorSpaceSubSpace.hpp"
#include "AbstractLinAlgPack/src/VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack/src/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/MultiVectorMutable.hpp"
#include "AbstractLinAlgPack/src/InnerProductDot.hpp"
#include "AbstractLinAlgPack/src/GenPermMatrixSlice.hpp"
#include "ThrowException.hpp"

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

VectorSpace::space_fcty_ptr_t
VectorSpace::small_vec_spc_fcty() const
{
	return MemMngPack::null;
}

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

VectorSpace::space_ptr_t
VectorSpace::space(
	const GenPermMatrixSlice  &P
	,BLAS_Cpp::Transp         P_trans
	) const
{
	const index_type
		dim = BLAS_Cpp::rows( P.rows(), P.cols(), P_trans );
	space_fcty_ptr_t  vec_spc_fcty = this->small_vec_spc_fcty();
	if(vec_spc_fcty.get())
		return vec_spc_fcty->create_vec_spc(dim);
	return MemMngPack::null;
}

// Overridden from AbstractFactory<>

VectorSpace::obj_ptr_t VectorSpace::create() const
{
	return this->create_member();
}

} // end namespace AbstractLinAlgPack
