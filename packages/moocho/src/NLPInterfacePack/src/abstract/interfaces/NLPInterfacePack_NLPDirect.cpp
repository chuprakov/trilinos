// /////////////////////////////////////////////////////
// NLPFirstOrderDirect.cpp
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

#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "Range1D.h"
#include "ThrowException.h"

namespace NLPInterfacePack {

// NLPFirstOrderDirect

size_type NLPFirstOrderDirect::r() const
{
	return this->con_decomp().size();
}

Range1D NLPFirstOrderDirect::var_dep() const
{
	return Range1D(1,m());
}
Range1D NLPFirstOrderDirect::var_indep() const
{
	return Range1D(m()+1,n());
}
Range1D NLPFirstOrderDirect::con_decomp() const
{
	return Range1D(1,m());
}

Range1D NLPFirstOrderDirect::con_undecomp() const
{
	return Range1D::Invalid;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_GcU() const
{
	return MemMngPack::null;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_Gh() const
{
	return MemMngPack::null;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_Uz() const
{
	return MemMngPack::null;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_Vz() const
{
	return MemMngPack::null;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_GcUD() const
{
	return MemMngPack::null;
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
NLPFirstOrderDirect::factory_GhD() const
{
	return MemMngPack::null;
}

NLP::vec_space_ptr_t NLPFirstOrderDirect::space_h() const
{
	return MemMngPack::null;
}

const VectorWithOp& NLPFirstOrderDirect::hl() const
{
	THROW_EXCEPTION( true, NoBounds, "NLPFirstOrderDirect::hl(), Error, default is for mI() == 0" );
	return xl(); // will never execute.
}

const VectorWithOp& NLPFirstOrderDirect::hu() const
{
	THROW_EXCEPTION( true, NoBounds, "NLPFirstOrderDirect::hl(), Error, default is for mI() == 0" );
	return xu(); // will never execute.
}

void NLPFirstOrderDirect::imp_calc_h(
	const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const
{
	assert(0); // Should never be called!
}

} // end namespace NLPIntefacePack
