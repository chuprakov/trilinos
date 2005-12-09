// //////////////////////////////////////////////////////////////////////
// VectorSpaceFactorySerial.cpp
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

#include "AbstractLinAlgPack_VectorSpaceFactorySerial.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpaceFactorySerial::VectorSpaceFactorySerial( const inner_prod_ptr_t& inner_prod )
	: VectorSpaceFactory(inner_prod)
{}

// Overridden from VectorSpaceFactory

VectorSpaceFactory::space_ptr_t
VectorSpaceFactorySerial::create_vec_spc(index_type dim) const
{
	return Teuchos::rcp(new VectorSpaceSerial(dim));
}

} // end namespace AbstractLinAlgPack
