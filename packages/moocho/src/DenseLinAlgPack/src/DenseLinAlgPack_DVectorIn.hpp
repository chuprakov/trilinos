// //////////////////////////////////////////////////////////////////////////////
// VectorIn.hpp
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

#ifndef VECTORIN_H
#define VECTORIN_H

#include "VectorInFunc.hpp"

namespace LinAlgPack {

///
/** Vector input stream operator.
  *
  * This operator function calls the function input(is,v,0).
  */
inline
std::istream& operator>>(std::istream& is, Vector& v)
{	return input(is,&v,(LinAlgPackIO::fmtflags)0); }

///
/** VectorSlice input stream operator.
  *
  * This operator function calls the function input(is,vs,0).
  */
inline
std::istream& operator>>(std::istream& is, VectorSlice& vs)
{	return input(is,&vs,(LinAlgPackIO::fmtflags)0); }

}	// end namespace LinAlgPack

#endif // VECTORIN_H
