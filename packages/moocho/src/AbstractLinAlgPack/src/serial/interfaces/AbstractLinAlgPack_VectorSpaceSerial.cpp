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

#include "SparseLinAlgPack/include/VectorSpaceSerial.h"
#include "SparseLinAlgPack/include/VectorWithOpMutableDense.h"
#include "SparseLinAlgPack/include/MultiVectorMutableDense.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/GenPermMatrixSlice.h"
#include "LinAlgPack/include/VectorClass.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

VectorSpaceSerial::VectorSpaceSerial( size_type dim )
{
	initialize(dim);
}

void VectorSpaceSerial::initialize( size_type dim )
{
	dim_ = dim;
}

bool VectorSpaceSerial::is_compatible(const VectorSpace& a_vec_space ) const
{
	return this->dim() == a_vec_space.dim();
}

// Overridden from VectorSpace

index_type VectorSpaceSerial::dim() const
{
	return dim_;
}

VectorSpace::space_ptr_t
VectorSpaceSerial::clone() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp( new VectorSpaceSerial( dim_	) );
}

VectorSpace::vec_mut_ptr_t
VectorSpaceSerial::create_member() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new VectorWithOpMutableDense(dim_));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceSerial::create_members(size_type num_vecs) const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new MultiVectorMutableDense(dim_,num_vecs));
}

VectorSpace::space_ptr_t
VectorSpaceSerial::sub_space(const Range1D& rng_in) const
{
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
	return MemMngPack::rcp( new VectorSpaceSerial( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ) ) ); 
}

} // end namespace SparseLinAlgPack
