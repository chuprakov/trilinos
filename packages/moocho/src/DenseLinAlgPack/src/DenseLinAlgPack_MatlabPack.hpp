// /////////////////////////////////////////////////////
// MatlabPack.h

#ifndef MATLAB_PACK_H
#define MATLAB_PACK_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {
namespace MatlabPack {

/** @name MatlabPack.
  *
  * This package contains functions that allow integration with 
  * Matlab.
  */
//@{

/** @name Output matrices and vectors in Matlab readable format.
  *
  * These function output vectors and matrices with enought digits
  * to reproduce the exact same floating point numbers.
  */
//@{

///
/** Output a VectorSlice.
  *
  * The VectorSlice is output in the following format:
  *
  \begin{verbatim}
	name = [ vs(1); vs(2); ... vs(n); ];
  \end{verbatim}
  *
  * Above #n = vs.size()# and #'# is appended to the end if #trans != BLAS_Cpp::no_trans#.
  * Also, a newline character #\n# is appended to the end after #']'#.
  */
std::ostream& out( std::ostream& o, const char* name, const VectorSlice& vs
	, BLAS_Cpp::Transp trans = BLAS_Cpp::no_trans );

///
/** Output a GenMatrixSlice.
  *
  * The GenMatrixSlice is output in the following format:
  *
  \begin{verbatim}
	name = [
	gms(1,1), gms(1,2), ... gms(1,n);
	gms(2,1), gms(2,2), ... gms(2,n);
	...       ...       ... ...
	gms(m,1), gms(m,2), ... gms(m,n);
	 ];
  \end{verbatim}
  *
  * Above #m = gms.rows()#, #n = gms.cols()# and #'# is appended to the end
  * if #trans != BLAS_Cpp::no_trans#.
  * Also, a newline character #\n# is appended to the end
  */
std::ostream& out( std::ostream& o, const char* name, const GenMatrixSlice& gms
	, BLAS_Cpp::Transp trans = BLAS_Cpp::no_trans );

//@}

//@}

}	// end namespace MatlabPack
}	// end namespace LinAlgPack

#endif	// MATLAB_PACK_H