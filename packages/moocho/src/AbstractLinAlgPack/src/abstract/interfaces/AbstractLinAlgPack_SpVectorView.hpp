// ////////////////////////////////////////////////////////////////////
// SpVectorView.h
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

#include "SpVectorClass.h"

namespace AbstractLinAlgPack {

///
/** Create an RTOp_SubVector view object from a SpVectorSlice object.
 */
RTOp_SubVector sub_vec_view(
	const SpVectorSlice&   sv
	,const Range1D&        rng = Range1D()
	);

} // end namespace AbstractLinAlgPack

#endif // SP_VECTOR_VIEW_H
