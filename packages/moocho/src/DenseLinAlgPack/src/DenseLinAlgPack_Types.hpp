// ///////////////////////////////////////////////////////////////////
// LinAlgPackTypes.h

#ifndef LINALGPACK_TYPES_H
#define LINALGPACK_TYPES_H

#include "BLAS_CppTypes.h"
#include "LinAlgPackOptions.h"

namespace LinAlgPack {

/** @name {\bf LinAlgPack Type Declarations}.
  *
  * These are forward declarations of the types used with the LinAlgPack
  * package (namespace).  In addition the BLAS_Cpp enumerations
  * \Ref{Transp}, \Ref{Side}, \Ref{Uplo}, and \Ref{Diag} and there values
  * are avalible using the qualifier #BLAS_Cpp#.
  */

//@{

///
using BLAS_Cpp::rows;
///
using BLAS_Cpp::cols;
///
using BLAS_Cpp::trans_not;

/// Enumeration for returning the amount of overlap between two objects
enum EOverLap { NO_OVERLAP = 0, SOME_OVERLAP, SAME_MEM };	

///
class IVector;
///
class Range1D;
///
class Vector;
///
class VectorSlice;
///
class TransVectorSlice;
///
class GenMatrix;
///
class GenMatrixSlice;
///
class TransGenMatrixSlice;
///
class tri_ele_gms;
///
class tri_gms;
///
class sym_gms;

//@}

}  // namespace LinAlgPack

#endif // LINALGPACK_TYPES_H