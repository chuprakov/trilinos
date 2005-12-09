// /////////////////////////////////////////////////////
// NLPDirect.cpp
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

#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_Range1D.hpp"
#include "Teuchos_TestForException.hpp"

namespace NLPInterfacePack {

// NLPDirect

void NLPDirect::set_factories(
	const mat_sym_fcty_ptr_t             &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	)
{
	factory_transDtD_ = factory_transDtD;
	factory_S_        = factory_S;
}

size_type NLPDirect::r() const
{
	return this->con_decomp().size();
}

Range1D NLPDirect::var_dep() const
{
	return Range1D(1,m());
}
Range1D NLPDirect::var_indep() const
{
	return Range1D(m()+1,n());
}
Range1D NLPDirect::con_decomp() const
{
	return Range1D(1,m());
}

Range1D NLPDirect::con_undecomp() const
{
	return Range1D::Invalid;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_GcU() const
{
	return Teuchos::null;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_Uz() const
{
	return Teuchos::null;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_GcUD() const
{
	return Teuchos::null;
}

const NLPDirect::mat_sym_fcty_ptr_t
NLPDirect::factory_transDtD() const
{
	return factory_transDtD_;
}
	
const NLPDirect::mat_sym_nonsing_fcty_ptr_t
NLPDirect::factory_S() const
{
	return factory_S_;
}

void NLPDirect::initialize(bool test_setup)
{
	NLPObjGrad::initialize(test_setup);
}

} // end namespace NLPIntefacePack
