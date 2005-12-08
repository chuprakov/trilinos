// //////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_InFormatDecl.hpp
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

#ifndef LINALGPACK_IN_FORMAT_DECL_H
#define LINALGPACK_IN_FORMAT_DECL_H

#include "DenseLinAlgPack_IOFormat.hpp"

namespace DenseLinAlgPack {

///
/* * Input stream operator function for bound_format objects.
  *
  * This template function performs the following tasks.
  * \begin{enumeration}
  * <li> Saves the formating state of #is#
  * <li> Sets the formating state of #is# to that stored in #bf.f()#
  * <li> Calls #input(is, bf.obj(), bf.extra_flags().flags())#
  * <li> Resets the streams formating state to its original.
  * \end{enumeration}
  *
  * The original formating state of #is# is preserved even if an exception
  * is thrown.
  */
template<class T>
std::istream& operator>>(std::istream& is, const LinAlgPackIO::bound_format<T>& bf);

}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_IN_FORMAT_DECL_H
