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
#include "LinAlgPack/include/VectorClass.h"

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
	namespace rcp = MemMngPack;
	return rcp::rcp( new VectorSpaceSerial( dim_	) );
}

VectorSpace::vec_mut_ptr_t
VectorSpaceSerial::create_member() const
{
	namespace rcp = MemMngPack;
	return rcp::rcp(new VectorWithOpMutableDense(dim_));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceSerial::create_members(size_type num_vecs) const
{
	namespace rcp = MemMngPack;
	return rcp::rcp(new MultiVectorMutableDense(dim_,num_vecs));
}

} // end namespace SparseLinAlgPack
