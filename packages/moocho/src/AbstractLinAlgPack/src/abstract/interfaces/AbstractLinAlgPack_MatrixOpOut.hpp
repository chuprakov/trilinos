// ///////////////////////////////////////////////////////////
// MatrixWithOpOut.h

#ifndef MATRIX_WITH_OP_OUT_H
#define MATRIX_WITH_OP_OUT_H

#include <iosfwd>

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

///
/** Output operator for \Ref{MatrixWithOp} objects.
 */
inline
std::ostream& operator<<( std::ostream& o, const MatrixWithOp& M )
{
	return M.output(o);
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_OUT_H
