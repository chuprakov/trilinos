// ///////////////////////////////////////////////////////////////////////////////////////
// BLAS_Cpp.h
//
// C++ overloads for BLAS kernals (element type removed from name and enum for operations)
//
// ToDo:  7/7/98:  Make more portable by using the types in FortranTypes.

#ifndef BLAS_CPP_OVERLOADS_DECLARATIONS_H
#define BLAS_CPP_OVERLOADS_DECLARATIONS_H

#include <stdexcept>

#include "BLAS_C_Decl.h"
#include "BLAS_CppTypes.h"

// Overloaded BLAS wrappers.
// The naming convention is the Fortran BLAS name minus the type prefix.
namespace BLAS_Cpp {

using namespace BLAS_C_Decl;

/** @name Option Arguments
 * These are emumerations that are used with the overloaded C++ BLAS declarations to replace the
 * error prone use of characters for specifying options
 * @memo enumerations (enum)
 */
//@{

///
const char SideChar[]	= {'L'	, 'R'			};
///
const char TransChar[]	= {'N'	, 'T'	, 'C'	};
///
const char UploChar[]	= {'U'	, 'L'			};
///
const char DiagChar[]	= {'U'	, 'N'			};

//@}

/** @name C++ BLAS Function Declarations
 * These are overloaded C++ functions that have removed the element type from the name
 * of the BLAS functions and use enumerations for the options arguments.
 */
//@{	

// ///////////////////////////////////////////////////////////////////////////////////////////
/** @name Level 1 BLAS (vector-vector operations) */
//@{	

/** @name Generate plane rotation */
//@{

///
inline
void rotg(double& a, double& b, double& c, double& s) {
	FORTRAN_FUNC_CALL(DROTG)(a,b,c,s);
}

//@}	 	

/** @name Generate modified plane rotation */
//@{

///
inline
void rotmg(double& d1, double& d2, double& a, const double& b, double* param) {
	FORTRAN_FUNC_CALL(DROTMG)(d1, d2, a, b, param);
}
//@}
 
/** @name Apply plane rotation */
//@{

///
inline
void rot(const int& N, double* X, const int& INCX, double* Y, const int& INCY
	, const double& C, const double& S)
{
	FORTRAN_FUNC_CALL(DROT)(N, X, INCX, Y, INCY, C, S);
}
//@}

/** @name  Apply modified plane rotation */
//@{

/// 
inline
void rot(const int& N, double* X, const int& INCX, double* Y, const int& INCY
	, const double* PARAM)
{
	FORTRAN_FUNC_CALL(DROTM)(N, X, INCX, Y, INCY, PARAM);
}

//@}

/** @name  Interchange vectors */
//@{

///
inline
void swap(const int& N, double* X, const int& INCX, double* Y, const int& INCY)
{
	FORTRAN_FUNC_CALL(DSWAP)(N, X, INCX, Y, INCY);
}		 	
//@}

/** @name  Vector scaling */
//@{

/// 
inline
void scal(const int& N, const double& ALPHA, double* X, const int& INCX)
{
	FORTRAN_FUNC_CALL(DSCAL)(N, ALPHA, X, INCX);
}
//@}

/** @name Vector copy */
//@{

/// 
inline
void copy(const int& N, const double* X, const int& INCX, double* Y, const int& INCY)
{
	FORTRAN_FUNC_CALL(DCOPY)(N, X, INCX, Y, INCY);
}
//@}

/** @name  y = a*x + y */
//@{

///
inline
void axpy(const int& N, const double& A, const double* X, const int& INCX, double* Y
	, const int& INCY)
{
	FORTRAN_FUNC_CALL(DAXPY)(N, A, X, INCX, Y, INCY);
}
//@}

/** @name  Dot product */
//@{

///
inline
double dot(const int& N, const double* X, const int& INCX, const double* Y, const int& INCY)
{
	return FORTRAN_FUNC_CALL(DDOT)(N, X, INCX, Y, INCY);
}

//@}

/** @name  2-Norm */
//@{

///
inline
double nrm2(const int& N, const double* X, const int& INCX)
{
	return FORTRAN_FUNC_CALL(DNRM2)(N, X, INCX);
}

//@}

/** @name  1-Norm */
//@{

///
inline
double asum(const int& N, const double* X, const int& INCX)
{
	return FORTRAN_FUNC_CALL(DASUM)(N, X, INCX);
}

//@}

/** @name  Inifinity-Norm */
//@{

///
inline
double iamax(const int& N, const double* X, const int& INCX)
{
	return FORTRAN_FUNC_CALL(IDAMAX)(N, X, INCX);
}

//@}

//		end Level-1 BLAS
//@}	

// /////////////////////////////////////////////////
/** @name Level-2 BLAS (matrix-vector operations) */
//@{	

/** @name General rectangular matrix-vector products */
//@{

///
inline
void gemv(Transp transa, int m, int n, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	FORTRAN_FUNC_CALL(DGEMV)(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}
			 
//@}	

/** @name General band matrix-vector products */
//@{

///
inline
void gbmv(Transp transa, int m, int n, int kl, int ku, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	FORTRAN_FUNC_CALL(DGBMV)(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
			 	
//@}

/** @name Hermitian matrix-vector products */
//@{

			 	
//@}

/** @name Hermitian band matrix-vector products */
//@{

//@}

/** @name Hermitian packed matrix-vector products */
//@{


//@}

/** @name Symmetric matrix-vector products */
//@{

///
inline
void symv(Uplo uplo, int n, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	FORTRAN_FUNC_CALL(DSYMV)(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
			
//@}

/** @name Symmetric band matrix-vector products */
//@{

///
inline
void sbmv(Uplo uplo, int n, int k, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	FORTRAN_FUNC_CALL(DSBMV)(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}

//@}

/** @name Symmetric packed matrix-vector products */
//@{

///
inline
void spmv(Uplo uplo, int n, double alpha, const double* pap
	, const double* x, int incx, double beta, double* py, int incy)
{
	FORTRAN_FUNC_CALL(DSPMV)(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}

//@}

/** @name Triangular matrix-vector products */
//@{

///
inline
void trmv(Uplo uplo, Transp trans, Diag diag, int n, const double* pa
	, int lda, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTRMV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}

//@}

/** @name Triangular band matrix-vector products */
//@{

///
inline
void tbmv(Uplo uplo, Transp trans, Diag diag, int n, int k, const double* pa
	, int lda, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTBMV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

//@}

/** @name Triangular packed matrix-vector products */
//@{

///
inline
void tpmv(Uplo uplo, Transp trans, Diag diag, int n, const double* pap
	, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTPMV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}

//@}

/** @name Triangular equation solve */
//@{

///
inline
void trsv(Uplo uplo, Transp trans, Diag diag, int n, const double* pa
	, int lda, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTRSV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}

//@}

/** @name Triangular band equation solve */
//@{

///
inline
void tbsv(Uplo uplo, Transp trans, Diag diag, int n, int k, const double* pa
	, int lda, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTBSV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

//@}

/** @name Triangular packed equation solve */
//@{

///
inline
void tpsv(Uplo uplo, Transp trans, Diag diag, int n, const double* pap
	, double* px, int incx)
{
	FORTRAN_FUNC_CALL(DTPSV)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}

//@}

/** @name General rank-1 update */
//@{

///
inline
void ger(int m, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pa, int lda)
{
	FORTRAN_FUNC_CALL(DGER)(m, n, alpha, px, incx, py, incy, pa, lda);
}

//@}

/** @name Hermitian rank-1 update */
//@{

//@}

/** @name Hermitian packed rank-1 update */
//@{

//@}

/** @name Hermitian rank-2 update */
//@{

//@}

/** @name Hermitian packed rank-2 update */
//@{

//@}

/** @name Symmetric rank-1 update */
//@{

///
inline
void syr(Uplo uplo, int n, double alpha, const double* px
	, int incx, double* pa, int lda)
{
	FORTRAN_FUNC_CALL(DSYR)(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}

//@}

/** @name Symmetric packed rank-1 update */
//@{

///
inline
void spr(Uplo uplo, int n, double alpha, const double* px
	, int incx, double* pap)
{
	FORTRAN_FUNC_CALL(DSPR)(&UploChar[uplo], n, alpha, px, incx, pap);
}

//@}

/** @name Symmetric rank-2 update */
//@{

///
inline
void syr2(Uplo uplo, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pa, int lda)
{
	FORTRAN_FUNC_CALL(DSYR2)(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}

//@}

/** @name Symmetric packed rank-2 update */
//@{

///
inline
void spr2(Uplo uplo, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pap)
{
	FORTRAN_FUNC_CALL(DSPR2)(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
}

//@}

//		end Level 2 BLAS
//@}
	
// /////////////////////////////////////////
/** @name Level 3 BLAS (matrix-matrix operations) */
//@{	

/** @name General rectangular matrix-matrix product */
//@{

///
inline
void gemm(Transp transa, Transp transb, int m, int n, int k, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	FORTRAN_FUNC_CALL(DGEMM)(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}

//@}

/** @name Symmetric matrix-matrix product */
//@{

///
inline
void symm(Side side, Uplo uplo, int m, int n, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	FORTRAN_FUNC_CALL(DSYMM)(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}

//@}

/** @name Hermitian matrix-matrix product */
//@{

//@}

/** @name Symmetric rank-k update */
//@{

///
inline
void syrk(Uplo uplo, Transp trans, int n, int k, double alpha, const double* pa
	, int lda, double beta, double* pc, int ldc)
{
	FORTRAN_FUNC_CALL(DSYRK)(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}

//@}

/** @name Hermitian rank-k update */
//@{

//@}

/** @name Symmetric rank-2k update */
//@{

///
inline
void syr2k(Uplo uplo, Transp trans, int n, int k, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	FORTRAN_FUNC_CALL(DSYR2K)(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}

//@}

/** @name Hermitian rank-2k update */
//@{

//@}

/** @name Triangular matrix-matrix product */
//@{

///
inline
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, double alpha
	, const double* pa, int lda, double* pb, int ldb)
{
	FORTRAN_FUNC_CALL(DTRMM)(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}

//@}

/** @name Solution of triangular system */
//@{

///
inline
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, double alpha
	, const double* pa, int lda, double* pb, int ldb)
{
	FORTRAN_FUNC_CALL(DTRSM)(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
		 	 	
//@}

//		end Level 3 BLAS
//@}	

//		end overloaded functions
//@}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_OVERLOADS_DECLARATIONS_H