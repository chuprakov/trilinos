// //////////////////////////////////////////////////////////////////
// AbstractFactory.hpp
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

#ifndef ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_H
#define ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_H

#include "ref_count_ptr.hpp"

namespace MemMngPack {

///
/** Simple, universal "Abstract Factory" interface for the dynamic
 * creation of objects.
 */
template<class T>
class AbstractFactory {
public:

#ifndef DOXYGEN_COMPILE
	///
	typedef MemMngPack::ref_count_ptr<T>   obj_ptr_t;
#endif

	///
	virtual ~AbstractFactory() {}

	///
	/** Create an object of type T returned as a smart reference
	 * counting pointer object.
	 */
	virtual obj_ptr_t create() const = 0;

}; // class AbstractFactory

} // end MemMngPack

#endif // ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_H
