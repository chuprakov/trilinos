// //////////////////////////////////////////////////////////////////////////////
// LinAlgPackInFormatDef.h
//
// Template definition file.

#ifndef LINALGPACK_IN_FORMAT_DEF_H
#define LINALGPACK_IN_FORMAT_DEF_H

#include "LinAlgPackInFormatDecl.h"

namespace LinAlgPack {

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

}	// end namespace LinAlgPack

#endif // LINALGPACK_IN_FORMAT_DEF_H
