// /////////////////////////////////////////////////////////////////////////////
// PermIn.h

#ifndef PERM_IN_H
#define PERM_IN_H

#include <istream>

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

///
/** Input stream operator for IVector used as a permutation array.
  *
  * The input format is:\\
  * #size#\\
  * #i1 i2 i3 ... isize#\\
  */
std::istream& operator>>(std::istream& istrm, IVector& perm);

}	// end namespace LinAlgPack

#endif // PERM_IN_H
