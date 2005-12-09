// //////////////////////////////////////////////////////////////
// BasisSystem.cpp
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

#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "Thyra_Range1D.hpp"

namespace AbstractLinAlgPack {

BasisSystem::BasisSystem(
	const mat_sym_fcty_ptr_t             &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	)
{
	this->initialize(factory_transDtD,factory_S);
}

void BasisSystem::initialize(
	const mat_sym_fcty_ptr_t             &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	)
{
	factory_transDtD_ = factory_transDtD;
	factory_S_        = factory_S;
}

Range1D BasisSystem::equ_decomp() const
{
	const size_type r = this->var_dep().size();
	return r ? Range1D(1,r) : Range1D::Invalid;
}

Range1D BasisSystem::equ_undecomp() const
{
	return Range1D::Invalid;
}

const BasisSystem::mat_fcty_ptr_t BasisSystem::factory_GcUP() const
{
	return Teuchos::null;
}

const BasisSystem::mat_sym_fcty_ptr_t
BasisSystem::factory_transDtD() const
{
	return factory_transDtD_;
}
	
const BasisSystem::mat_sym_nonsing_fcty_ptr_t
BasisSystem::factory_S() const
{
	return factory_S_;
}

} // end namespace AbstractLinAlgPack
