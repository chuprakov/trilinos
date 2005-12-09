// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SpVectorView.hpp
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

#ifndef SP_VECTOR_VIEW_H
#define SP_VECTOR_VIEW_H

#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "RTOpPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Create an <tt>RTOpPack::SparseSubVector</tt> view object from a
 * <tt>SpVectorSlice</tt> object.
 */
RTOpPack::SparseSubVector sub_vec_view(
	const SpVectorSlice&   sv
	,const Range1D&        rng = Range1D()
	);

} // end namespace AbstractLinAlgPack

#endif // SP_VECTOR_VIEW_H
