// ////////////////////////////////////////////////////////////////////////////
// VectorSpaceSerial.cpp
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

#include "SparseLinAlgPack/src/VectorSpaceSerial.hpp"
#include "SparseLinAlgPack/src/VectorSpaceFactorySerial.hpp"
#include "SparseLinAlgPack/src/VectorWithOpMutableDense.hpp"
#include "SparseLinAlgPack/src/MultiVectorMutableDense.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/GenPermMatrixSlice.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "ThrowException.hpp"

#ifdef _DEBUG
#define CLASS_MEMBER_PTRS \
const VectorSpaceSerial  *_this = this; \
const size_type *_dim = &dim_;
#else
#define CLASS_MEMBER_PTRS
#endif

namespace SparseLinAlgPack {

VectorSpaceSerial::VectorSpaceSerial( size_type dim )
{
	CLASS_MEMBER_PTRS
	initialize(dim);
}

void VectorSpaceSerial::initialize( size_type dim )
{
	CLASS_MEMBER_PTRS
	dim_ = dim;
}

bool VectorSpaceSerial::is_compatible(const VectorSpace& a_vec_space ) const
{
	CLASS_MEMBER_PTRS
	return this->dim() == a_vec_space.dim();
}

// Overridden from VectorSpace

index_type VectorSpaceSerial::dim() const
{
	CLASS_MEMBER_PTRS
	return dim_;
}

VectorSpace::space_fcty_ptr_t
VectorSpaceSerial::small_vec_spc_fcty() const
{
	CLASS_MEMBER_PTRS
	return MemMngPack::rcp(new VectorSpaceFactorySerial());
}

VectorSpace::space_ptr_t
VectorSpaceSerial::clone() const
{
	CLASS_MEMBER_PTRS
	namespace mmp = MemMngPack;
	return mmp::rcp( new VectorSpaceSerial( dim_	) );
}

VectorSpace::vec_mut_ptr_t
VectorSpaceSerial::create_member() const
{
	CLASS_MEMBER_PTRS
	namespace mmp = MemMngPack;
	return mmp::rcp(new VectorWithOpMutableDense(dim_));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceSerial::create_members(size_type num_vecs) const
{
	CLASS_MEMBER_PTRS
	namespace mmp = MemMngPack;
	return mmp::rcp(new MultiVectorMutableDense(dim_,num_vecs));
}

VectorSpace::space_ptr_t
VectorSpaceSerial::sub_space(const Range1D& rng_in) const
{
	CLASS_MEMBER_PTRS
	namespace mmp = MemMngPack;
	const size_type this_dim = this->dim();
	const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"VectorSpaceSerial::sub_view(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim );
#endif
	if( rng == Range1D(1,this_dim) )
		return mmp::rcp( this, false );
	return mmp::rcp( new VectorSpaceSerial( rng.size() ) ); 
}

VectorSpace::space_ptr_t
VectorSpaceSerial::space(
	const GenPermMatrixSlice  &P
	,BLAS_Cpp::Transp         P_trans
	) const
{
	CLASS_MEMBER_PTRS
	return MemMngPack::rcp( new VectorSpaceSerial( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ) ) ); 
}

} // end namespace SparseLinAlgPack
