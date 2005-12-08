// //////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_InFormatDef.hpp
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
//
// Template definition file.

#ifndef LINALGPACK_IN_FORMAT_DEF_H
#define LINALGPACK_IN_FORMAT_DEF_H

#include <istream>

#include "DenseLinAlgPack_InFormatDecl.hpp"

namespace DenseLinAlgPack {

template<class T>
std::istream& operator>>(std::istream& is, const LinAlgPackIO::bound_format<T>& bf) {
	using LinAlgPackIO::ios_format_memento;

	ios_format_memento old_format = ios_format_memento::save_format(is);

	try {
		bf.f().set_format(is);
		input( is, &const_cast< LinAlgPackIO::bound_format<T>&>(bf).obj()
			   , bf.f().extra_flags().flags() );
	}
	catch(...) {
		old_format.set_format(is);
		throw;
	}

	old_format.set_format(is);
	return is;
}

}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_IN_FORMAT_DEF_H
