// //////////////////////////////////////////////////////////////////////////////
// VectorIn.h

#ifndef VECTORIN_H
#define VECTORIN_H

#include "VectorInFunc.h"

namespace LinAlgPack {

///
/** Vector input stream operator.
  *
  * This operator function calls the function input(is,v,0).
  */
inline
std::istream& operator>>(std::istream& is, Vector& v)
{	return input(is,&v,(LinAlgPackIO::fmtflags)0); }

///
/** VectorSlice input stream operator.
  *
  * This operator function calls the function input(is,vs,0).
  */
inline
std::istream& operator>>(std::istream& is, VectorSlice& vs)
{	return input(is,&vs,(LinAlgPackIO::fmtflags)0); }

}	// end namespace LinAlgPack

#endif // VECTORIN_H