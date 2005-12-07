// ///////////////////////////////////////////////////////////////////////////
// BLAS_Cpp_Types.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//
// Types for BLAS_Cpp namespace
//

#ifndef BLAS_CPP_TYPES_H
#define BLAS_CPP_TYPES_H

#include "Moocho_ConfigDefs.hpp"

namespace BLAS_Cpp {

/** \defgroup BLAS_Cpp_grp Basic C++/Fortran BLAS declarations/utilities.
 * \ingroup Misc_grp
 *
 * These are declarations that are useful in BLAS (Basic Linear Algebra
 * Subroutines) types of codes.  They provide some foundation for C++/Fortran
 * compatibility.
 */
//@{

/// Size type
typedef size_t size_type;

/** \defgroup BLAS_Cpp_grp_enums BLAS Enumerations */
//@{

/// SIDE
enum Side{
	left
	,right
};
/// TRANS
enum Transp{
	no_trans     ///< Not transposed
	,trans       ///< Transposed
	,conj_trans  ///< Conjugate transpose
};
/// UPLO
enum Uplo {upper, lower};
/// Return the opposite of Uplo argument
inline Uplo operator!(Uplo uplo)
{	return uplo == upper ? lower : upper; }
/// DIAG
enum Diag {unit, nonunit};

//@}

/** \defgroup BLAS_Cpp_grp_helper Helper functions */
//@{

/// Return Transp given a bool
inline Transp bool_to_trans(bool return_trans) {
	return return_trans ? trans : no_trans;
}

/// Returns true if _trans == trans
inline bool trans_to_bool(Transp _trans) {
	return ( _trans == trans );
}

/// Return the opposite of the transpose argument
inline Transp operator!(Transp _trans)	// does not work with conj_trans
{	return _trans == no_trans ? trans : no_trans; }

/// Return the opposite of the transpose argument
inline Transp trans_not(Transp _trans)	// does not work with conj_trans
{	return _trans == no_trans ? trans : no_trans; }

/// Return the transpose of the transpose argument
inline Transp trans_trans(Transp _trans1, Transp _trans2) // does not work with conj_trans
{	return _trans1 == _trans2 ? no_trans : trans; }

/// Give a string name to Transp value
inline const char*  trans_to_string(Transp _trans)
{
	return _trans == no_trans ? "no_trans" : "trans";
}

/// Return rows of a possible transposed matrix
inline size_type rows(size_type rows, size_type cols
	, BLAS_Cpp::Transp _trans)
{
	return _trans == BLAS_Cpp::no_trans ? rows : cols;
}

/// Return columns of a possible transposed matrix
inline size_type cols(size_type rows, size_type cols
	, BLAS_Cpp::Transp _trans)
{
	return _trans == BLAS_Cpp::no_trans ? cols : rows;
}

//@}

//@}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_TYPES_H
