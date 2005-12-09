// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_COOMatrixOut.hpp
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

#ifndef COO_MATRIX_OUT_H
#define COO_MATRIX_OUT_H

#include "AbstractLinAlgPack_COOMatrixOutFunc.hpp"

namespace AbstractLinAlgPack {

///
/** Output stream operator for COOMatrix.
  *
  * Calls the function:
  *
  * std::ostream& output(std::ostream& o, const COOMatrix& coom)
  */
std::ostream& operator<<(std::ostream& o, const COOMatrix& coom);

// /////////////////////////////////////////////////////////////
// Inline definition

inline std::ostream& operator<<(std::ostream& o, const COOMatrix& coom) {
	return output(o,coom);
}

} // end namespace AbstractLinAlgPack

#endif // COO_MATRIX_OUT_H
