// ///////////////////////////////////////////////////////////////
// MultiVectorMutable.cpp
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

#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"

namespace AbstractLinAlgPack {

MultiVectorMutable::vec_ptr_t MultiVectorMutable::row(index_type i) const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(
		const_cast<MultiVectorMutable*>(this)->row(i));
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::col(index_type j) const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(
		const_cast<MultiVectorMutable*>(this)->col(j));
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::diag(int k) const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(
		const_cast<MultiVectorMutable*>(this)->diag(k));
}

} // end namespace AbstractLinAlgPack
