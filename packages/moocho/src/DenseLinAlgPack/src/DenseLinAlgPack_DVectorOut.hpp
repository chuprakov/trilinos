// //////////////////////////////////////////////////////////////////////////////////////
// VectorOut.h
//
// Output stream operator for Vector

#ifndef VECTOROUT_H
#define VECTOROUT_H

#include "VectorOutFunc.h"

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