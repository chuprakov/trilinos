// ///////////////////////////////////////////////////////////////
// InnerProductDot.cpp
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

#include "AbstractLinAlgPack/src/InnerProductDot.hpp"
#include "AbstractLinAlgPack/src/VectorWithOp.hpp"
#include "RTOpStdOpsLib/src/RTOp_ROp_dot_prod.h"
#include "RTOpPack/src/RTOpCppC.hpp"

namespace AbstractLinAlgPack {

value_type InnerProductDot::inner_prod(
	const VectorWithOp& v1, const VectorWithOp& v2
	) const
{
	RTOpPack::RTOpC         op;
	RTOpPack::ReductTarget  reduct_obj;
	RTOp_ROp_dot_prod_construct( &op.op() );
	op.reduct_obj_create(&reduct_obj);
	const VectorWithOp* vecs[1] = { &v2 };
	v1.apply_reduction(op,1,vecs,0,NULL,reduct_obj.obj());
	RTOp_value_type val = RTOp_ROp_dot_prod_val(reduct_obj.obj());
	return val;
}

} // end namespace AbstractLinAlgPack
