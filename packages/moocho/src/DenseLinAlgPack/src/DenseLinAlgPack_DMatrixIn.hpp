// ///////////////////////////////////////////////////////////////////////////////////////
// GenMatrixIn.h

#ifndef GENMATRIX_IN_H
#define GENMATRIX_IN_H

#include "GenMatrixInFunc.h"

namespace LinAlgPack {

///
/** GenMatrix input stream operator.
  *
  * This operator function calls the function input(is,gm,0).
  */
inline
std::istream& operator>>(std::istream& is, GenMatrix& gm)
{	return input(is,&gm,0); }

///
/** GenMatrixSlice input stream operator.
  *
  * This operator function calls the function input(is,gms,0).
  */
inline
std::istream& operator>>(std::istream& is, GenMatrixSlice& gms)
{	return input(is,&gms,0); }

}	// end namespace LinAlgPack

#endif // GENMATRIX_IN_H