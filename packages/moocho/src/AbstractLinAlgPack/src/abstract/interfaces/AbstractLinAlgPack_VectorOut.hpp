// ///////////////////////////////////////////////////////////
// VectorWithOpOut.h

#ifndef VECTOR_WITH_OP_OUT_H
#define VECTOR_WITH_OP_OUT_H

#include <iosfwd>

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

///
/** Output operator for \Ref{VectorWithOp} objects.
 */
inline
std::ostream& operator<<( std::ostream& o, const VectorWithOp& v )
{
	return v.output(o);
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_OUT_H
