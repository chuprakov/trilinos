// ///////////////////////////////////////////////////////////////////////////////////////
// GenMatrixOut.h

#ifndef GENMATRIX_OUT_H
#define GENMATRIX_OUT_H

#include "GenMatrixOutFunc.h"

namespace LinAlgPack {

///
/** GenMatrixSlice output stream operator.
  *
  * This operator function calls the function output(os,gms,0).
  */
inline std::ostream& operator<<(std::ostream& os, const GenMatrixSlice& gms) {
	return output(os, gms, 0);
}

}	// end namespace LinAlgPack

// ////////////////////////////////////
// Inline function definitions

//inline std::ostream& LinAlgPack::operator<<(std::ostream& os, const GenMatrixSlice& gms) {
//	return output(os, gms, 0);
//}

#endif // GENMATRIX_OUT_H