// //////////////////////////////////////////////////////////////////////////////////////
// VectorOut.hpp
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
// Output stream operator for Vector

#ifndef VECTOROUT_H
#define VECTOROUT_H

#include "VectorOutFunc.hpp"

namespace LinAlgPack {

///
/** VectorSlice output stream operator.
  *
  * This operator function calls the function output(os,vs,0).
  */
inline std::ostream& operator<<(std::ostream& os, const VectorSlice& vs) {
	return output(os, vs, (LinAlgPackIO::fmtflags)(0));
}

}	// end namespace LinAlgPack

// ////////////////////////////////////
// Inline function definitions

//inline std::ostream& LinAlgPack::operator<<(std::ostream& os, const VectorSlice& vs) {
//	return output(os, vs, 0);
//}

#endif // VECTOROUT_H
