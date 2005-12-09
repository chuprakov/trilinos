// //////////////////////////////////////////////////////////////////////
// VectorSpaceFactory.cpp
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

#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack_InnerProduct.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpaceFactory::~VectorSpaceFactory()
{}

VectorSpaceFactory::VectorSpaceFactory( const inner_prod_ptr_t& inner_prod )
{
	this->inner_prod(inner_prod);
}

void VectorSpaceFactory::inner_prod( const inner_prod_ptr_t& inner_prod )
{
	inner_prod_ = inner_prod;
}

const VectorSpaceFactory::inner_prod_ptr_t
VectorSpaceFactory::inner_prod() const
{
	return inner_prod_;
}

} // end namespace AbstractLinAlgPack
