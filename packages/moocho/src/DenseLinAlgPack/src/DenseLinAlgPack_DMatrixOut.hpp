// ///////////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DMatrixOut.hpp
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

#ifndef GENMATRIX_OUT_H
#define GENMATRIX_OUT_H

#include "DenseLinAlgPack_DMatrixOutFunc.hpp"

namespace DenseLinAlgPack {

///
/* * DMatrixSlice output stream operator.
  *
  * This operator function calls the function output(os,gms,0).
  */
inline std::ostream& operator<<(std::ostream& os, const DMatrixSlice& gms) {
	return output(os, gms, (LinAlgPackIO::fmtflags)(0));
}

}	// end namespace DenseLinAlgPack

// ////////////////////////////////////
// Inline function definitions

//inline std::ostream& DenseLinAlgPack::operator<<(std::ostream& os, const DMatrixSlice& gms) {
//	return output(os, gms, 0);
//}

#endif // GENMATRIX_OUT_H
