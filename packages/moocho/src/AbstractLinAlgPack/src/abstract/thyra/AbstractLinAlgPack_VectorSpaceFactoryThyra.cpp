// ///////////////////////////////////////////////////////////////
// VectorSpaceFactoryThyra.cpp
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

#include "VectorSpaceFactoryThyra.hpp"
#include "VectorSpaceThyra.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

VectorSpaceFactoryThyra::VectorSpaceFactoryThyra()
{}

VectorSpaceFactoryThyra::VectorSpaceFactoryThyra(
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> >& thyra_vec_spc_fcty
	)
{
	this->initialize(thyra_vec_spc_fcty);
}

void VectorSpaceFactoryThyra::initialize(
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> >& thyra_vec_spc_fcty
	)
{
	TEST_FOR_EXCEPTION(
		thyra_vec_spc_fcty.get()==NULL, std::invalid_argument
		,"VectorSpaceFactoryThyra::initialize(thyra_vec_spc_fcty): Error!"
		);
	thyra_vec_spc_fcty_ = thyra_vec_spc_fcty;
}

Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> > 
VectorSpaceFactoryThyra::set_uninitialized()
{
	Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> > tmp_thyra_vec_spc_fcty = thyra_vec_spc_fcty_;
	thyra_vec_spc_fcty_ = Teuchos::null;
	return tmp_thyra_vec_spc_fcty;
}

// Overridden from VectorSpaceFactory

VectorSpaceFactory::space_ptr_t
VectorSpaceFactoryThyra::create_vec_spc(index_type dim) const
{
	return Teuchos::rcp(new VectorSpaceThyra(thyra_vec_spc_fcty_->createVecSpc(dim)));
}

} // end namespace AbstractLinAlgPack
