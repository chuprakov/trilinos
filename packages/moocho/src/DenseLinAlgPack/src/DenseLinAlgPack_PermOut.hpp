// /////////////////////////////////////////////////////////////////////////////
// PermOut.h

#ifndef PERM_OUT_H
#define PERM_OUT_H

#include <ostream>

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

///
/** Output stream operator for IVector used as a permutation array.
  *
  * The output format is:\\
  * #size#\\
  * #i1 i2 i3 ... isize#\\
  */
std::ostream& operator<<(std::ostream& o, const IVector& perm);

}	// end namespace LinAlgPack

#endif // PERM_OUT_H