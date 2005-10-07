// ///////////////////////////////////////////////////////////////
// InnerProductThyra.cpp
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

#include <stdexcept>

#include "InnerProductThyra.hpp"
#include "VectorMutableThyra.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

InnerProductThyra::InnerProductThyra()
{}

InnerProductThyra::InnerProductThyra(
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
	)
{
	this->initialize(thyra_vec_spc);
}

void InnerProductThyra::initialize(
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
	)
{
	TEST_FOR_EXCEPTION(
		thyra_vec_spc.get()==NULL, std::invalid_argument
		,"InnerProductThyra::initialize(thyra_vec_spc): Error!"
		);
	thyra_vec_spc_ = thyra_vec_spc;
}

Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > 
InnerProductThyra::set_uninitialized()
{
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > tmp_thyra_vec_spc = thyra_vec_spc_;
	thyra_vec_spc_ = Teuchos::null;
	return tmp_thyra_vec_spc;
}

// Overridden from InnerProduct

value_type InnerProductThyra::inner_prod(const Vector& v1, const Vector& v2) const
{
	using Teuchos::dyn_cast;
	return thyra_vec_spc_->scalarProd(
		*dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
		,*dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
		);
}

} // end namespace AbstractLinAlgPack
