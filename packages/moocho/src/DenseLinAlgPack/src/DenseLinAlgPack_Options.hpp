// ////////////////////////////////////////////////////////////////////////////////////////////////
// linalg_options.h
//
// Options for LinAlgPack compilation
//

#ifndef LINALG_OPTIONS_H
#define LINALG_OPTIONS_H

#include "LinAlgPackDebugAcronyms.h"
#include "Misc/include/fortran_types.h"

/** @name {\bf LinAlgPack Options}.
  *
  * The header file LinAlgPackOptions.h contains the defines for several macros that
  * determine how the library is built.  The user should comment out any
  * macros that her or she does not want to be defined.  The definition of
  * these macros cause the library code to assert the preconditions documented
  * for each of the member and non-member functions and throw the listed exceptions
  * if they are not satisfied.  Precondtions are supposed to be the 
  * responcibility of the client code so the user may only want to define
  * these macros during debugging for better program verification.
  * If the user checks all of the preconditions listed in this documentation for the calls
  * to all functions then the checks performed by the library are redundant.
  */
//@{

///
/** If defined the library code checks to see if subscripts are in bounds for element access
  * an subregion indexing.  If the preconditions for the subscripting operations are
  * not satisfied then the listed exceptions will be thrown.
  */
#define LINALGPACK_CHECK_RANGE

///
/** If defined the library code checks to see if the sizes of rhs arguments in expressions are compatible.
  * The exception std::length_error will be thrown if rhs sizes are not compatible.
  */
#define LINALGPACK_CHECK_RHS_SIZES

///
/** If defined the library code checks to see if VectorSlice and GenMatrixSlice objects have valid constructions.
  * If they do not have valid constructions then an exception will be thrown.  The operation of these
  * checks may depend on the definition of the macro \Ref{LINALGPACK_CHECK_RANGE}.
  */
#define LINALGPACK_CHECK_SLICE_SETUP

#if defined LINALGPACK_NO_CHECKS

// Turn them all off.
#undef LINALGPACK_CHECK_RANGE
#undef LINALGPACK_CHECK_RHS_SIZES
#undef LINALGPACK_CHECK_SLICE_SETUP

#endif

namespace LinAlgPack{
/// Typedef for the value type of elements that is used for the library.
typedef FortranTypes::f_dbl_prec		value_type;
/// Typedef for the indice type of elements that are used by the library
typedef FortranTypes::f_int				indice_type;
/// Typedef for the size type of elements that are used by the library
typedef	size_t							size_type;

}

//@}

#endif