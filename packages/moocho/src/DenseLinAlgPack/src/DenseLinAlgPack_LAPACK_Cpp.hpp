// ///////////////////////////////////////////////////////
// DenseLinAlgPack_LAPACK_Cpp.hpp
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

#ifndef LAPACK_CPP_H
#define LAPACK_CPP_H

#include "DenseLinAlgPack_LAPACK_C_Decl.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"

// Cpp Declarations for calling LAPACK functions that
// use function overloading to remove the floating point
// typeds from the names and use enumerations to replace
// the char arguments.

namespace LAPACK_Cpp {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;
typedef FortranTypes::f_logical		f_logical;

// xPOTRF

inline
///
void potrf(  BLAS_Cpp::Uplo uplo
	, const f_int& n, f_dbl_prec* A, const f_int& lda
	, f_int* info )
{
	LAPACK_C_Decl::dpotrf( BLAS_Cpp::UploChar[uplo]
		,n ,A ,lda ,info );
}

// xGEQRF

inline
///
void geqrf( const f_int& m
	, const f_int& n, f_dbl_prec* A, const f_int& lda
	, f_dbl_prec* tau, f_dbl_prec* work
	, const f_int& lwork, f_int* info  )
{
	LAPACK_C_Decl::dgeqrf(m,n,A,lda,tau,work,lwork,info);
}

// xORMQR

inline
///
void ormqr( BLAS_Cpp::Side side, BLAS_Cpp::Transp trans
	, const f_int& m, const f_int& n
	, const f_int& k, const f_dbl_prec* A, const f_int& lda
	, const f_dbl_prec* tau, f_dbl_prec* C, const f_int& ldc
	, f_dbl_prec* work, const f_int& lwork, f_int* info )
{
	LAPACK_C_Decl::dormqr( BLAS_Cpp::SideChar[side]
		, BLAS_Cpp::TransChar[trans], m, n, k, A, lda
		, tau, C, ldc, work, lwork, info );
	
}

// xSYTRF

inline
//
void sytrf( BLAS_Cpp::Uplo uplo, const f_int& n, f_dbl_prec A[]
	, const f_int& lda, f_int ipiv[], f_dbl_prec work[], const f_int& lwork
	, f_int* info )
{
	LAPACK_C_Decl::dsytrf( BLAS_Cpp::UploChar[uplo]
		, n, A, lda, ipiv, work, lwork, info );
}

// xSYTRS

inline
///
void sytrs( BLAS_Cpp::Uplo uplo
	, const f_int& n, const f_int& nrhs, const f_dbl_prec A[]
	, const f_int& lda, const f_int ipiv[], f_dbl_prec B[]
	, const f_int& ldb, f_int* info )
{
	LAPACK_C_Decl::dsytrs( BLAS_Cpp::UploChar[uplo]
		, n, nrhs, A, lda, ipiv, B, ldb, info );
}

// xGETRF

inline
//
void getrf(
	const f_int& m, const f_int& n, f_dbl_prec A[]
	,const f_int& lda, f_int ipiv[], f_int* info )
{
	LAPACK_C_Decl::dgetrf( m, n, A, lda, ipiv, info );
}

// xGETRS

inline
///
void getrs(
	BLAS_Cpp::Transp trans
	,const f_int& n, const f_int& nrhs, const f_dbl_prec A[]
	, const f_int& lda, const f_int ipiv[], f_dbl_prec B[]
	, const f_int& ldb, f_int* info )
{
	LAPACK_C_Decl::dgetrs(
		BLAS_Cpp::TransChar[trans], n, nrhs, A, lda, ipiv, B, ldb, info
		);
}

}	// end namespace LAPACK_Cpp

#endif	// LAPACK_CPP_H
