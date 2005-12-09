// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_COOMatrixIn.hpp
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

#ifndef COO_MATRIX_IN_H
#define COO_MATRIX_IN_H

#include <istream>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Inputstream stream operator for COOMatrix.
  *
  * The format for the imput is:
  *
  * #m  n  nz#\\
  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
  *
  * In the above format, each non-zero element is given as a three item pair:
  * value of the non-zero element, row indice (1-based) of the non-zero element,
  * and the column indice (1-based) of the non-zero element.  There must
  * be no spaces between the numbers and the \':\' charachter and there must
  * be at least one whitespace character between elements.
  */
std::istream& operator>>(std::istream& istrm, COOMatrix& coom);

// Inline definition
inline std::istream& operator>>(std::istream& istrm, COOMatrix& coom) {
	coom.initialize(istrm); return istrm;
}

} // end namespace AbstractLinAlgPack

#endif // COO_MATRIX_IN_H
