// //////////////////////////////////////////////////////////////////////////////
// COOMatrixOut.h

#ifndef COO_MATRIX_OUT_H
#define COO_MATRIX_OUT_H

#include "COOMatrixOutFunc.h"

namespace SparseLinAlgPack {

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

} // end namespace SparseLinAlgPack

#endif // COO_MATRIX_OUT_H