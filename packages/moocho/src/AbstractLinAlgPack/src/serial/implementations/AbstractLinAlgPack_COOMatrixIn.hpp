// //////////////////////////////////////////////////////////////////////////////
// COOMatrixIn.h

#ifndef COO_MATRIX_IN_H
#define COO_MATRIX_IN_H

#include <istream>

#include "SparseLinAlgPackTypes.h"

namespace SparseLinAlgPack {

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

} // end namespace SparseLinAlgPack

#endif // COO_MATRIX_IN_H
