// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorSpaceFactorySerial.hpp
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

#ifndef VECTOR_SPACE_FACTORY_SERIAL_H
#define VECTOR_SPACE_FACTORY_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * ToDo: Finish documentation!
 */
class VectorSpaceFactorySerial
	: public AbstractLinAlgPack::VectorSpaceFactory // doxygen needs full name
{
public:

	///
	VectorSpaceFactorySerial( const inner_prod_ptr_t& inner_prod = Teuchos::null );

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	space_ptr_t create_vec_spc(index_type dim) const;

	//@}
	
}; // end class VectorSpaceFactorySerial

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_SPACE_FACTORY_SERIAL_H
