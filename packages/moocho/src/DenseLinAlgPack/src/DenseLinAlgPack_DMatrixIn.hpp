// ///////////////////////////////////////////////////////////////////////////////////////
// GenMatrixIn.h
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

#ifndef GENMATRIX_IN_H
#define GENMATRIX_IN_H

#include "GenMatrixInFunc.h"

namespace LinAlgPack {

///
/** GenMatrix input stream operator.
  *
  * This operator function calls the function input(is,gm,0).
  */
inline
std::istream& operator>>(std::istream& is, GenMatrix& gm)
{	return input(is,&gm,(LinAlgPackIO::fmtflags)0); }

///
/** GenMatrixSlice input stream operator.
  *
  * This operator function calls the function input(is,gms,0).
  */
inline
std::istream& operator>>(std::istream& is, GenMatrixSlice& gms)
{	return input(is,&gms,(LinAlgPackIO::fmtflags)0); }

}	// end namespace LinAlgPack

#endif // GENMATRIX_IN_H
