#ifndef TSFLAPACK_H
#define TSFLAPACK_H

#include "TSFConfig.h"

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
	void dgemv_(const char* transFlag, const int* nRows, const int* nCols, 
							const double* alpha, const double* A, const int* lda,
							const double* x, const int* incx, const double* beta,
							double* y, int* incx);

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
		};

	/** \ingroup Utilities
	 * Specialization of TSFBlas to double 
	 */

	template<> class TSFBlas<double>
		{
		public:
			inline static void copy(const int* n, const double* x, const int* incx, 
															double* y, const int* incy)
				{dcopy_(n, x, incx, y, incy);}
			
			
			inline static void axpy(const int* n, const double* a, const double* x, 
															int* incx, double* y, int* incy)
				{daxpy_(n, a, x, incx, y, incy);}
				
			inline static double dot(const int* n, const double* x, const int* incx, 
															 const double* y, const int* incy)
				{return ddot_(n, x, incx, y, incy);}
			
			inline static double nrm2(const int* n, const double* x, const int* incx)
				{return dnrm2_(n, x, incx);}
			
			inline	static double scal(const int* n, const double* a, 
																 const double* x, const int* incx)
				{return dscal_(n, a, x, incx);}
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

			inline	static float scal(const int* n, const float* a, 
																const float* x, const int* incx)
				{return sscal_(n, a, x, incx);}
		};

}


#endif /* end of Blas.h */
