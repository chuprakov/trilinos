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

#include "MultiVectorMutable.hpp"
#include "VectorMutable.hpp"
#include "WorkspacePack.hpp"

namespace AbstractLinAlgPack {

// Sub-view methods

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
	assert(0); // ToDo: return a MultiVectorMutableSubView object.
	// Note that the MultiVectorMutableSubView class should derive from
	// MultiVectorSubView.
	return MemMngPack::null;
}

// Overriddend form MultiVector

MultiVectorMutable::vec_ptr_t MultiVectorMutable::row(index_type i) const
{
	return const_cast<MultiVectorMutable*>(this)->row(i);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::col(index_type j) const
{
	return const_cast<MultiVectorMutable*>(this)->col(j);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::diag(int k) const
{
	return const_cast<MultiVectorMutable*>(this)->diag(k);
}

MultiVectorMutable::multi_vec_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	return const_cast<MultiVectorMutable*>(this)->mv_sub_view(row_rng,col_rng);
}

} // end namespace AbstractLinAlgPack
