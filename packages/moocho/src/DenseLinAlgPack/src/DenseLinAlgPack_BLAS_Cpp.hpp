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
void rotg(float& a, float& b, float& c, float& s) {
	srotg(a,b,c,s);
}		 	
///
inline
void rotg(double& a, double& b, double& c, double& s) {
	drotg(a,b,c,s);
}

//@}	 	

/** @name Generate modified plane rotation */
//@{

///
inline
void rotmg(float& d1, float& d2, float& a, const float& b, float* param) {
	srotmg(d1, d2, a, b, param);
}
///
inline
void rotmg(double& d1, double& d2, double& a, const double& b, double* param) {
	drotmg(d1, d2, a, b, param);
}
//@}
 
/** @name Apply plane rotation */
//@{

///
inline
void rot(const int& N, float* X, const int& INCX, float* Y, const int& INCY
	, const float& C, const float& S)
{
	srot(N, X, INCX, Y, INCY, C, S);
}
///
inline
void rot(const int& N, double* X, const int& INCX, double* Y, const int& INCY
	, const double& C, const double& S)
{
	drot(N, X, INCX, Y, INCY, C, S);
}
//@}

/** @name  Apply modified plane rotation */
//@{

///
inline
void rot(const int& N, float* X, const int& INCX, float* Y, const int& INCY
	, const float* PARAM)
{
	srotm(N, X, INCX, Y, INCY, PARAM);
}
/// 
inline
void rot(const int& N, double* X, const int& INCX, double* Y, const int& INCY
	, const double* PARAM)
{
	drotm(N, X, INCX, Y, INCY, PARAM);
}

//@}

/** @name  Interchange vectors */
//@{

///
inline
void swap(const int& N, float* X, const int& INCX, float* Y, const int& INCY)
{
	sswap(N, X, INCX, Y, INCY);
}
///
inline
void swap(const int& N, double* X, const int& INCX, double* Y, const int& INCY)
{
	dswap(N, X, INCX, Y, INCY);
}		 	
///
inline
void swap(const int& N, complex<float>* X, const int& INCX, complex<float>* Y
	, const int& INCY)
{
	cswap(N, X, INCX, Y, INCY);
}		 	
///
inline
void swap(const int& N, complex<double>* X, const int& INCX, complex<double>* Y
	, const int& INCY)
{
	zswap(N, X, INCX, Y, INCY);
}
//@}

/** @name  Vector scaling */
//@{

/// 
inline
void scal(const int& N, const float& ALPHA, float* X, const int& INCX)
{
	sscal(N, ALPHA, X, INCX);
}
/// 
inline
void scal(const int& N, const double& ALPHA, double* X, const int& INCX)
{
	dscal(N, ALPHA, X, INCX);
}
/// 
inline
void scal(const int& N, const complex<float>& ALPHA, complex<float>* X, const int& INCX)
{
	cscal(N, ALPHA, X, INCX);
}
///
inline
void scal(const int& N, const complex<double>& ALPHA, complex<double>* X, const int& INCX)
{
	zscal(N, ALPHA, X, INCX);
}
/// 
inline
void scal(const int& N, const float& ALPHA, complex<float>* X, const int& INCX)
{
	csscal(N, ALPHA, X, INCX);
}
///
inline
void scal(const int& N, const double& ALPHA, complex<double>* X, const int& INCX)
{
	zdscal(N, ALPHA, X, INCX);
}
//@}

/** @name Vector copy */
//@{

/// 
inline
void copy(const int& N, const float* X, const int& INCX, float* Y, const int& INCY)
{
	scopy(N, X, INCX, Y, INCY);
}
/// 
inline
void copy(const int& N, const double* X, const int& INCX, double* Y, const int& INCY)
{
	dcopy(N, X, INCX, Y, INCY);
}
///
inline
void copy(const int& N, const complex<float>* X, const int& INCX, complex<float>* Y
	, const int& INCY)
{
	ccopy(N, X, INCX, Y, INCY);
}
///
inline
void copy(const int& N, const complex<double>* X, const int& INCX, complex<double>* Y
	, const int& INCY)
{
	zcopy(N, X, INCX, Y, INCY);
}
//@}

/** @name  y = a*x + y */
//@{

///
inline
void axpy(const int& N, const float& A, const float* X, const int& INCX, float* Y
	, const int& INCY)
{
	saxpy(N, A, X, INCX, Y, INCY);
}
///
inline
void axpy(const int& N, const double& A, const double* X, const int& INCX, double* Y
	, const int& INCY)
{
	daxpy(N, A, X, INCX, Y, INCY);
}
///
inline
void axpy(const int& N, const complex<float>& A, const complex<float>* X, const int& INCX
	, complex<float>* Y, const int& INCY) {caxpy(N, A, X, INCX, Y, INCY);
}
///
inline
void axpy(const int& N, const complex<double>& A, const complex<double>* X, const int& INCX
	, complex<double>* Y, const int& INCY) {zaxpy(N, A, X, INCX, Y, INCY);
}
//@}

/** @name  Dot product */
//@{

///
inline
float dot(const int& N, const float* X, const int& INCX, const float* Y, const int& INCY)
{
	return sdot(N, X, INCX, Y, INCY);
}
///
inline
double dot(const int& N, const double* X, const int& INCX, const double* Y, const int& INCY)
{
	return ddot(N, X, INCX, Y, INCY);
}
///
inline
double sdot(const int& N, const float* X, const int& INCX, const float* Y, const int& INCY)
{
	return dsdot(N, X, INCX, Y, INCY);
}
///
inline
complex<float> dotu(const int& N, const complex<float>* X, const int& INCX
	, const complex<float>* Y, const int& INCY) 
{
	complex<float> temp;
	cdotu(temp, N, X, INCX, Y, INCY);
	return temp;
}
///
inline
complex<double> dotu(const int& N, const complex<double>* X, const int& INCX
, const complex<double>* Y, const int& INCY)
{
	complex<double> temp;
	zdotu(temp, N, X, INCX, Y, INCY);
	return temp;
}
///
inline
complex<float> dotc(const int& N, const complex<float>* X, const int& INCX
	, const complex<float>* Y, const int& INCY)
{
	complex<float> temp;
	cdotc(temp, N, X, INCX, Y, INCY);
	return temp;
}
///
inline
complex<double> dotc(const int& N, const complex<double>* X, const int& INCX
	, const complex<double>* Y, const int& INCY)
{
	complex<double> temp;
	zdotc(temp, N, X, INCX, Y, INCY);
	return temp;
}
///
inline
float dot(const int& N, const float& B, const float* X, const int& INCX, const float* Y
	, const int& INCY)
{
	return sdsdot(N, B, X, INCX, Y, INCY);
}
//@}

/** @name  2-Norm */
//@{

///
inline
float nrm2(const int& N, const float* X, const int& INCX)
{
	return snrm2(N, X, INCX);
}
///
inline
double nrm2(const int& N, const double* X, const int& INCX)
{
	return dnrm2(N, X, INCX);
}
///
inline
complex<float> nrm2(const int& N, const complex<float>* X, const int& INCX)
{
	return scnrm2(N, X, INCX);
}
///
inline
complex<double> nrm2(const int& N, const complex<double>* X, const int& INCX)
{
	return dznrm2(N, X, INCX);
}

//@}

/** @name  1-Norm */
//@{

///
inline
float asum(const int& N, const float* X, const int& INCX)
{
	return sasum(N, X, INCX);
}
///
inline
double asum(const int& N, const double* X, const int& INCX)
{
	return dasum(N, X, INCX);
}
///
inline
complex<float> asum(const int& N, const complex<float>* X, const int& INCX)
{
	return scasum(N, X, INCX);
}
///
inline
complex<double> asum(const int& N, const complex<double>* X, const int& INCX)
{
	return dzasum(N, X, INCX);
}

//@}

/** @name  Inifinity-Norm */
//@{

///
inline
float iamax(const int& N, const float* X, const int& INCX)
{
	return isamax(N, X, INCX);
}
///
inline
double iamax(const int& N, const double* X, const int& INCX)
{
	return idamax(N, X, INCX);
}
///
inline
complex<float> iamax(const int& N, const complex<float>* X, const int& INCX)
{
	return icamax(N, X, INCX);
}
///
inline
complex<double>
iamax(const int& N, const complex<double>* X, const int& INCX)
{
	return izamax(N, X, INCX);
}

//@}

//		end Level 1 BLAS
//@}	

// /////////////////////////////////////////////////
/** @name Level 2 BLAS (matrix-vector operations) */
//@{	

/** @name General rectangular matrix-vector products */
//@{

///
inline
void gemv(Transp transa, int m, int n, float alpha, const float* pa
	, int lda, const float* x, int incx, float beta, float* py, int incy)
{
	sgemv(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gemv(Transp transa, int m, int n, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	dgemv(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gemv(Transp transa, int m, int n, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* x, int incx
	, const complex<float>& beta, complex<float>* py, int incy)
{
	cgemv(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gemv(Transp transa, int m, int n, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* x, int incx
	, const complex<double>& beta, complex<double>* py, int incy)
{
	zgemv(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}
			 
//@}	

/** @name General band matrix-vector products */
//@{

///
inline
void gbmv(Transp transa, int m, int n, int kl, int ku, float alpha, const float* pa
	, int lda, const float* x, int incx, float beta, float* py, int incy)
{
	sgbmv(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gbmv(Transp transa, int m, int n, int kl, int ku, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	dgbmv(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gbmv(Transp transa, int m, int n, int kl, int ku, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* x, int incx
	, const complex<float>& beta, complex<float>* py, int incy)
{
	cgbmv(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void gbmv(Transp transa, int m, int n, int kl, int ku, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* x, int incx
	, const complex<double>& beta, complex<double>* py, int incy)
{
	zgbmv(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
			 	
//@}

/** @name Hermitian matrix-vector products */
//@{

///
inline
void hemv(Uplo uplo, int n, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* x, int incx
	, const complex<float>& beta, complex<float>* py, int incy)
{
	chemv(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void hemv(Uplo uplo, int n, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* x, int incx
	, const complex<double>& beta, complex<double>* py, int incy)
{
	zhemv(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
				 	
//@}

/** @name Hermitian band matrix-vector products */
//@{

///
inline
void hbmv(Uplo uplo, int n, int k, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* x, int incx
	, const complex<float>& beta, complex<float>* py, int incy)
{
	chbmv(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void hbmv(Uplo uplo, int n, int k, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* x, int incx
	, const complex<double>& beta, complex<double>* py, int incy)
{
	zhbmv(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}

//@}

/** @name Hermitian packed matrix-vector products */
//@{

///
inline
void hpmv(Uplo uplo, int n, const complex<float>& alpha
	, const complex<float>* pap, const complex<float>* x, int incx
	, const complex<float>& beta, complex<float>* py, int incy)
{
	chpmv(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}		 	
///
inline
void hpmv(Uplo uplo, int n, const complex<double>& alpha
	, const complex<double>* pap, const complex<double>* x, int incx
	, const complex<double>& beta, complex<double>* py, int incy)
{
	zhpmv(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}
	
//@}

/** @name Symmetric matrix-vector products */
//@{

///
inline
void symv(Uplo uplo, int n, float alpha, const float* pa
	, int lda, const float* x, int incx, float beta, float* py, int incy)
{
	ssymv(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void symv(Uplo uplo, int n, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	dsymv(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
			
//@}

/** @name Symmetric band matrix-vector products */
//@{

///
inline
void sbmv(Uplo uplo, int n, int k, float alpha, const float* pa
	, int lda, const float* x, int incx, float beta, float* py, int incy)
{
	ssbmv(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}
///
inline
void sbmv(Uplo uplo, int n, int k, double alpha, const double* pa
	, int lda, const double* x, int incx, double beta, double* py, int incy)
{
	dsbmv(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}
//@}

/** @name Symmetric packed matrix-vector products */
//@{

///
inline
void spmv(Uplo uplo, int n, float alpha, const float* pap
	, const float* x, int incx, float beta, float* py, int incy)
{
	sspmv(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}
///
inline
void spmv(Uplo uplo, int n, double alpha, const double* pap
	, const double* x, int incx, double beta, double* py, int incy)
{
	dspmv(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}

//@}

/** @name Triangular matrix-vector products */
//@{

///
inline
void trmv(Uplo uplo, Transp trans, Diag diag, int n, const float* pa
	, int lda, float* px, int incx)
{
	strmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trmv(Uplo uplo, Transp trans, Diag diag, int n, const double* pa
	, int lda, double* px, int incx)
{
	dtrmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trmv(Uplo uplo, Transp trans, Diag diag, int n, const complex<float>* pa
	, int lda, complex<float>* px, int incx)
{
	ctrmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trmv(Uplo uplo, Transp trans, Diag diag, int n, const complex<double>* pa
	, int lda, complex<double>* px, int incx)
{
	ztrmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
//@}

/** @name Triangular band matrix-vector products */
//@{

///
inline
void tbmv(Uplo uplo, Transp trans, Diag diag, int n, int k, const float* pa
	, int lda, float* px, int incx)
{
	stbmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}
///
inline
void tbmv(Uplo uplo, Transp trans, Diag diag, int n, int k, const double* pa
	, int lda, double* px, int incx)
{
	dtbmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}
///
inline
void tbmv(Uplo uplo, Transp trans, Diag diag, int n, int k, const complex<float>* pa
	, int lda, complex<float>* px, int incx)
{
	ctbmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}
///
inline
void tbmv(Uplo uplo, Transp trans, Diag diag, int n, int k, const complex<double>* pa
	, int lda, complex<double>* px, int incx)
{
	ztbmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}
//@}

/** @name Triangular packed matrix-vector products */
//@{

///
inline
void tpmv(Uplo uplo, Transp trans, Diag diag, int n, const float* pap
	, float* px, int incx)
{
	stpmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpmv(Uplo uplo, Transp trans, Diag diag, int n, const double* pap
	, double* px, int incx)
{
	dtpmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpmv(Uplo uplo, Transp trans, Diag diag, int n, const complex<float>* pap
	, complex<float>* px, int incx)
{
	ctpmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpmv(Uplo uplo, Transp trans, Diag diag, int n, const complex<double>* pap
	, complex<double>* px, int incx)
{
	ztpmv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
//@}

/** @name Triangular equation solve */
//@{

///
inline
void trsv(Uplo uplo, Transp trans, Diag diag, int n, const float* pa
	, int lda, float* px, int incx)
{
	strsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trsv(Uplo uplo, Transp trans, Diag diag, int n, const double* pa
	, int lda, double* px, int incx)
{
	dtrsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trsv(Uplo uplo, Transp trans, Diag diag, int n, const complex<float>* pa
	, int lda, complex<float>* px, int incx)
{
	ctrsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
///
inline
void trsv(Uplo uplo, Transp trans, Diag diag, int n, const complex<double>* pa
	, int lda, complex<double>* px, int incx)
{
	ztrsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}
//@}

/** @name Triangular band equation solve */
//@{

///
inline
void tbsv(Uplo uplo, Transp trans, Diag diag, int n, int k, const float* pa
	, int lda, float* px, int incx)
{
	stbsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

///
inline
void tbsv(Uplo uplo, Transp trans, Diag diag, int n, int k, const double* pa
	, int lda, double* px, int incx)
{
	dtbsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

///
inline
void tbsv(Uplo uplo, Transp trans, Diag diag, int n, int k, const complex<float>* pa
	, int lda, complex<float>* px, int incx)
{
	ctbsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

///
inline
void tbsv(Uplo uplo, Transp trans, Diag diag, int n, int k, const complex<double>* pa
	, int lda, complex<double>* px, int incx)
{
	ztbsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}
//@}

/** @name Triangular packed equation solve */
//@{

///
inline
void tpsv(Uplo uplo, Transp trans, Diag diag, int n, const float* pap
	, float* px, int incx)
{
	stpsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpsv(Uplo uplo, Transp trans, Diag diag, int n, const double* pap
	, double* px, int incx)
{
	dtpsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpsv(Uplo uplo, Transp trans, Diag diag, int n, const complex<float>* pap
	, complex<float>* px, int incx)
{
	ctpsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
///
inline
void tpsv(Uplo uplo, Transp trans, Diag diag, int n, const complex<double>* pap
	, complex<double>* px, int incx)
{
	ztpsv(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}
//@}

/** @name General rank-1 update */
//@{

///
inline
void ger(int m, int n, float alpha, const float* px
	, int incx, const float* py, int incy, float* pa, int lda)
{
	sger(m, n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void ger(int m, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pa, int lda)
{
	dger(m, n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void geru(int m, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, const complex<float>* py, int incy, complex<float>* pa, int lda)
{
	cgeru(m, n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void geru(int m, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, const complex<double>* py, int incy, complex<double>* pa, int lda)
{
	zgeru(m, n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void gerc(int m, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, const complex<float>* py, int incy, complex<float>* pa, int lda)
{
	cgerc(m, n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void gerc(int m, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, const complex<double>* py, int incy, complex<double>* pa, int lda)
{
	zgerc(m, n, alpha, px, incx, py, incy, pa, lda);
}
//@}

/** @name Hermitian rank-1 update */
//@{

///
inline
void her(Uplo uplo, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, complex<float>* pa, int lda)
{
	cher(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}
///
inline
void her(Uplo uplo, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, complex<double>* pa, int lda)
{
	zher(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}
//@}

/** @name Hermitian packed rank-1 update */
//@{

///
inline
void hpr(Uplo uplo, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, complex<float>* pap)
{
	chpr(&UploChar[uplo], n, alpha, px, incx, pap);
}
///
inline
void hpr(Uplo uplo, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, complex<double>* pap)
{
	zhpr(&UploChar[uplo], n, alpha, px, incx, pap);
}
//@}

/** @name Hermitian rank-2 update */
//@{

///
inline
void her2(Uplo uplo, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, const complex<float>* py, int incy, complex<float>* pa, int lda)
{
	cher2(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}
///
inline
void her2(Uplo uplo, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, const complex<double>* py, int incy, complex<double>* pa, int lda)
{
	zher2(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}
//@}

/** @name Hermitian packed rank-2 update */
//@{

///
inline
void hpr2(Uplo uplo, int n, const complex<float>& alpha, const complex<float>* px
	, int incx, const complex<float>* py, int incy, complex<float>* pap)
{
	chpr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
}
///
inline
void hpr2(Uplo uplo, int n, const complex<double>& alpha, const complex<double>* px
	, int incx, const complex<double>* py, int incy, complex<double>* pap)
{
	zhpr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
}
//@}

/** @name Symmetric rank-1 update */
//@{

///
inline
void syr(Uplo uplo, int n, float alpha, const float* px
	, int incx, float* pa, int lda)
{
	ssyr(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}
///
inline
void syr(Uplo uplo, int n, double alpha, const double* px
	, int incx, double* pa, int lda)
{
	dsyr(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}
//@}

/** @name Symmetric packed rank-1 update */
//@{

///
inline
void spr(Uplo uplo, int n, float alpha, const float* px
	, int incx, float* pap)
{
	sspr(&UploChar[uplo], n, alpha, px, incx, pap);
}		 	
///
inline
void spr(Uplo uplo, int n, double alpha, const double* px
	, int incx, double* pap)
{
	dspr(&UploChar[uplo], n, alpha, px, incx, pap);
}
//@}

/** @name Symmetric rank-2 update */
//@{

///
inline
void syr2(Uplo uplo, int n, float alpha, const float* px
	, int incx, const float* py, int incy, float* pa, int lda)
{
	ssyr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}		 	
///
inline
void syr2(Uplo uplo, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pa, int lda)
{
	dsyr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}
//@}

/** @name Symmetric packed rank-2 update */
//@{

///
inline
void spr2(Uplo uplo, int n, float alpha, const float* px
	, int incx, const float* py, int incy, float* pap)
{
	sspr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
}
///
inline
void spr2(Uplo uplo, int n, double alpha, const double* px
	, int incx, const double* py, int incy, double* pap)
{
	dspr2(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
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
void gemm(Transp transa, Transp transb, int m, int n, int k, float alpha, const float* pa
	, int lda, const float* pb, int ldb, float beta, float* pc, int ldc)
{
	sgemm(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}
///
inline
void gemm(Transp transa, Transp transb, int m, int n, int k, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	dgemm(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}
///
inline
void gemm(Transp transa, Transp transb, int m, int n, int k, complex<float> alpha, const complex<float>* pa
	, int lda, const complex<float>* pb, int ldb, complex<float> beta, complex<float>* pc, int ldc)
{
	cgemm(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}		 	
///
inline
void gemm(Transp transa, Transp transb, int m, int n, int k, complex<double> alpha, const complex<double>* pa
	, int lda, const complex<double>* pb, int ldb, complex<double> beta, complex<double>* pc, int ldc)
{
	zgemm(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}
//@}

/** @name Symmetric matrix-matrix product */
//@{

///
inline
void symm(Side side, Uplo uplo, int m, int n, float alpha, const float* pa
	, int lda, const float* pb, int ldb, float beta, float* pc, int ldc)
{
	ssymm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void symm(Side side, Uplo uplo, int m, int n, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	dsymm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void symm(Side side, Uplo uplo, int m, int n, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* pb, int ldb
	, const complex<float>& beta, complex<float>* pc, int ldc)
{
	csymm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void symm(Side side, Uplo uplo, int m, int n, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* pb, int ldb
	, const complex<double>& beta, complex<double>* pc, int ldc)
{
	zsymm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
//@}

/** @name Hermitian matrix-matrix product */
//@{

///
inline
void hemm(Side side, Uplo uplo, int m, int n, const complex<float>& alpha
	, const complex<float>* pa, int lda, const complex<float>* pb, int ldb
	, const complex<float>& beta, complex<float>* pc, int ldc)
{
	chemm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void hemm(Side side, Uplo uplo, int m, int n, const complex<double>& alpha
	, const complex<double>* pa, int lda, const complex<double>* pb, int ldb
	, const complex<double>& beta, complex<double>* pc, int ldc)
{
	zhemm(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
//@}

/** @name Symmetric rank-k update */
//@{

///
inline
void syrk(Uplo uplo, Transp trans, int n, int k, float alpha, const float* pa
	, int lda, float beta, float* pc, int ldc)
{
	ssyrk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}
///
inline
void syrk(Uplo uplo, Transp trans, int n, int k, double alpha, const double* pa
	, int lda, double beta, double* pc, int ldc)
{
	dsyrk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}
///
inline
void syrk(Uplo uplo, Transp trans, int n, int k, complex<float> alpha
	, const complex<float>* pa, int lda, complex<float> beta, complex<float>* pc, int ldc)
{
	csyrk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}		 	
///
inline
void syrk(Uplo uplo, Transp trans, int n, int k, complex<double> alpha
	, const complex<double>* pa, int lda, complex<double> beta, complex<double>* pc, int ldc)
{
	zsyrk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}
//@}

/** @name Hermitian rank-k update */
//@{

///
inline
void herk(Uplo uplo, Transp trans, int n, int k, float alpha
	, const complex<float>* pa, int lda, float beta, complex<float>* pc, int ldc)
{
	cherk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}		 	
///
inline
void herk(Uplo uplo, Transp trans, int n, int k, double alpha
	, const complex<double>* pa, int lda, double beta, complex<double>* pc, int ldc)
{
	zherk(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}
//@}

/** @name Symmetric rank-2k update */
//@{

///
inline
void syr2k(Uplo uplo, Transp trans, int n, int k, float alpha, const float* pa
	, int lda, const float* pb, int ldb, float beta, float* pc, int ldc)
{
	ssyr2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void syr2k(Uplo uplo, Transp trans, int n, int k, double alpha, const double* pa
	, int lda, const double* pb, int ldb, double beta, double* pc, int ldc)
{
	dsyr2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void syr2k(Uplo uplo, Transp trans, int n, int k, complex<float> alpha
	, const complex<float>* pa, int lda, const complex<float>* pb, int ldb
	, complex<float> beta, complex<float>* pc, int ldc)
{
	csyr2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void syr2k(Uplo uplo, Transp trans, int n, int k, complex<double> alpha
	, const complex<double>* pa, int lda, const complex<double>* pb, int ldb
	, complex<double> beta, complex<double>* pc, int ldc)
{
	zsyr2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
//@}

/** @name Hermitian rank-2k update */
//@{

///
inline
void her2k(Uplo uplo, Transp trans, int n, int k, complex<float> alpha
	, const complex<float>* pa, int lda, const complex<float>* pb, int ldb
	, float beta, complex<float>* pc, int ldc)
{
	cher2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
///
inline
void her2k(Uplo uplo, Transp trans, int n, int k, complex<double> alpha
	, const complex<double>* pa, int lda, const complex<double>* pb, int ldb
	, double beta, complex<double>* pc, int ldc)
{
	zher2k(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}
//@}

/** @name Triangular matrix-matrix product */
//@{

///
inline
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, float alpha
	, const float* pa, int lda, float* pb, int ldb)
{
	strmm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, double alpha
	, const double* pa, int lda, double* pb, int ldb)
{
	dtrmm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n
	, complex<float> alpha, const complex<float>* pa, int lda
	, complex<float>* pb, int ldb)
{
	ctrmm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n
	, complex<double> alpha, const complex<double>* pa, int lda
	, complex<double>* pb, int ldb)
{
	ztrmm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
//@}

/** @name Solution of triangular system */
//@{

///
inline
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, float alpha
	, const float* pa, int lda, float* pb, int ldb)
{
	strsm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n, double alpha
	, const double* pa, int lda, double* pb, int ldb)
{
	dtrsm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n
	, complex<float> alpha, const complex<float>* pa, int lda
	, complex<float>* pb, int ldb)
{
	ctrsm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
///
inline
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, int m, int n
	, complex<double> alpha, const complex<double>* pa, int lda
	, complex<double>* pb, int ldb)
{
	ztrsm(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}
		 	 	
//@}

//		end Level 3 BLAS
//@}	

//		end overloaded functions
//@}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_OVERLOADS_DECLARATIONS_H