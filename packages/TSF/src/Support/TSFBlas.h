#ifndef TSFBLAS_H
#define TSFBLAS_H

#include "TSFDefs.h"
#include "TSFTimeMonitor.h"

extern "C"
{
  /* double precision methods */
  void dcopy_(const int* n, const double* x, const int* incx,
              double* y, const int* incy);
  void daxpy_(const int* n, const double* a, const double* x,
              int* incx, double* y, int* incy);
  double ddot_(const int* n, const double* x, const int* incx,
               const double* y, const int* incy);
  double dnrm2_(const int* n, const double* x, const int* incx);
  double dscal_(const int* n, const double* a,
                const double* x, const int* incx);
  /* matrix-vector multiply: y = alpha*A*x + beta*y */
  void dgemv_(const char* transFlag, const int* nRows, const int* nCols,
              const double* alpha, const double* A, const int* lda,
              const double* x, const int* incx, const double* beta,
              double* y, int* incy);

  /* factor a general matrix */
  void dgetrf_(const int* nRows, const int* nCols, double* A,
               const int* lda, int* iPiv, int* info);

  /* backsolve a factored matrix */
  void dgetrs_(const char* transFlag, const int* nEqn, const int* nRHS,
               const double* A, const int* lda, const int* ipiv,
               double* b, const int* ldb, int* info);

  /* single precision methods */
  void scopy_(const int* n, const float* x, const int* incx,
              float* y, const int* incy);
  void saxpy_(const int* n, const float* a, const float* x,
              int* incx, float* y, int* incy);
  float sdot_(const int* n, const float* x, const int* incx,
              const float* y, const int* incy);
  float snrm2_(const int* n, const float* x, const int* incx);
  float sscal_(const int* n, const float* a,
               const float* x, const int* incx);
  void sgemv_(const char* transFlag, const int* nRows, const int* nCols,
              const float* alpha, const float* A, const int* lda,
              const float* x, const int* incx, const float* beta,
              float* y, int* incy);
  void sgetrf_(const int* nRows, const int* nCols, float* A,
               const int* lda, int* iPiv, int* info);
  /* backsolve a factored matrix */
  void sgetrs_(const char* transFlag, const int* nEqn, const int* nRHS,
               const float* A, const int* lda, const int* ipiv,
               float* b, const int* ldb, int* info);
}

namespace TSF
{
  using std::string;

  /** \ingroup Utilities
   * Wrapper for selected BLAS functions, templated on TSF real type.
   */
  template<class Real> class TSFBlas
    {
    public:
      /** copy vector x into vector y */
      static void copy(const int* n, const Real* x, const int* incx,
                       Real* y, const int* incy);

      /** axpy: y += a*x, returning y by reference argument.
       * Argument a is a scalar, x and y are vectors */
      static void axpy(const int* n, const Real* a, const Real* x,
                       int* incx, Real* y, int* incy);

      /** return dot product of vectors x and y */
      static Real dot(const int* n, const Real* x, const int* incx,
                      const Real* y, const int* incy);

      /** return 2-norm of vector x */
      static Real nrm2(const int* n, const Real* x, const int* incx);

      /** multiply vector x by a scalar constant a */
      static Real scal(const int* n, const Real* a,
                       const Real* x, const int* incx);

      /** do matrix-vector multiply: y = alpha*A*x + beta*y */
      static void gemv(const char* transFlag, const int* nRows,
                       const int* nCols,
                       const Real* alpha, const Real* A, const int* lda,
                       const Real* x, const int* incx, const Real* beta,
                       Real* y, int* incx);
      /* factor a general matrix */
      static void getrf(const int* nRows, const int* nCols, Real* A,
                        const int* lda, int* iPiv, int* info);

      /* backsolve a factored matrix */
      static void getrs(const char* transFlag, const int* nEqn,
                        const int* nRHS,
                        const Real* A, const int* lda, const int* ipiv,
                        Real* b, const int* ldb, int* info);



    };

  /** \ingroup Utilities
   * Specialization of TSFBlas to double
   */

  template<> class TSFBlas<double>
    {
    public:
      inline static void copy(const int* n, const double* x, const int* incx,
                              double* y, const int* incy)
        {TSFTimeMonitor t(blasTimer()); dcopy_(n, x, incx, y, incy);}


      inline static void axpy(const int* n, const double* a, const double* x,
                              int* incx, double* y, int* incy)
        {TSFTimeMonitor t(blasTimer()); daxpy_(n, a, x, incx, y, incy);}

      inline static double dot(const int* n, const double* x, const int* incx,
                               const double* y, const int* incy)
        {TSFTimeMonitor t(blasTimer()); return ddot_(n, x, incx, y, incy);}

      inline static double nrm2(const int* n, const double* x, const int* incx)
        {TSFTimeMonitor t(blasTimer()); return dnrm2_(n, x, incx);}

      inline  static double scal(const int* n, const double* a,
                                 const double* x, const int* incx)
        {TSFTimeMonitor t(blasTimer()); return dscal_(n, a, x, incx);}

      /** do matrix-vector multiply: y = alpha*A*x + beta*y */
      inline static void gemv(const char* transFlag, const int* nRows,
                              const int* nCols,
                              const double* alpha, const double* A,
                              const int* lda,
                              const double* x, const int* incx,
                              const double* beta,
                              double* y, int* incy)
        {dgemv_(transFlag, nRows, nCols, alpha, A, lda, x, incx,
                beta, y, incy);}

      /* factor a general matrix */
      inline static void getrf(const int* nRows, const int* nCols, double* A,
                               const int* lda, int* iPiv, int* info)
        {dgetrf_(nRows, nCols, A, lda, iPiv, info);}


      /* backsolve a factored matrix */
      inline static void getrs(const char* transFlag, const int* nEqn,
                               const int* nRHS,
                               const double* A, const int* lda,
                               const int* ipiv,
                               double* b, const int* ldb, int* info)
        {dgetrs_(transFlag, nEqn, nRHS, A, lda, ipiv, b, ldb, info);}

      static TSFTimer& blasTimer()
        {
          static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("BLAS calls");
          return *timer;
        }
    };


  /** \ingroup Utilities
   * Specialization of TSFBlas to float
   */
  template<> class TSFBlas<float>
    {
    public:
      inline static void copy(const int* n, const float* x, const int* incx,
                              float* y, const int* incy)
        {scopy_(n, x, incx, y, incy);}


      inline static void axpy(const int* n, const float* a, const float* x,
                              int* incx, float* y, int* incy)
        {saxpy_(n, a, x, incx, y, incy);}

      inline static float dot(const int* n, const float* x, const int* incx,
                              const float* y, const int* incy)
        {return sdot_(n, x, incx, y, incy);}

      inline static float nrm2(const int* n, const float* x, const int* incx)
        {return snrm2_(n, x, incx);}

      inline  static float scal(const int* n, const float* a,
                                const float* x, const int* incx)
        {return sscal_(n, a, x, incx);}
    };

}


#endif /* end of Blas.h */
