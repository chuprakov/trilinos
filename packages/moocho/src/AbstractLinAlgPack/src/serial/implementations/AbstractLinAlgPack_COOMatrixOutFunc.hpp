// //////////////////////////////////////////////////////////////////////////////
// COOMatrixOutFunc.h

#ifndef COO_MATRIX_OUT_FUNC_H
#define COO_MATRIX_OUT_FUNC_H

#include <ostream>

#include "SparseLinAlgPackTypes.h"

namespace SparseLinAlgPack {

///
/** Output stream function for COOMatrix.
  *
  * The format for the output is:
  *
  * #m  n  nz#\\
  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
  *
  * In the above format, each non-zero element is given as a three item pair:
  * value of the non-zero element, row indice (1-based) of the non-zero element,
  * and the column indice (1-based) of the non-zero element.  There are
  * no spaces between the numbers and the \':\' charachter and there is
  * one whitespace character between elements.  There is a new line character
  * added after all of the elements have been output.
  */
std::ostream& output(std::ostream& o, const COOMatrix& coom);

} // end namespace SparseLinAlgPack

#endif // COO_MATRIX_OUT_FUNC_H
