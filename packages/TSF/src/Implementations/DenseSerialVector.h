#ifndef DENSESERIALVECTOR_H
#define DENSESERIALVECTOR_H

#include "TSFConfig.h"
#include "TSFError.h"
#include "TSFTimer.h"
#include <stdlib.h>
#include <iostream.h>
#include <string.h>



namespace TSF
{
	using std::string;
	

	/**\ingroup DenseSerial
	 * Serial vector with math operations
	 */

	class DenseSerialVector
		{
		public:
			/** {\bf Constructors, Destructors, and Assignment Operator} */
			//@{
			/** Empty ctor */
			DenseSerialVector() : x_(0), n_(0) {;}
			/** Create a vector of length n */
			inline DenseSerialVector(int n);
			/** Create a vector of length n, and fill it with value */
			inline DenseSerialVector(int n, TSFReal value);
			/** Create a vector of length n, assuming responsibility for a C array */
			inline DenseSerialVector(int n, const TSFReal* cArray);
			/** Copy ctor */
			inline DenseSerialVector(const DenseSerialVector& other);
			/**  dtor */
			inline  ~DenseSerialVector();
			/** assignment op */
			DenseSerialVector& operator=(const DenseSerialVector& other);
			//@}

			/**  {\bf element access with optional bounds checking} */
			//@{
			/** read/write access */
			inline TSFReal& operator[](int i);
			/** read-only access */
			inline const TSFReal& operator[](int i) const ;

			//@}
	
			/** {\bf sizing } */
			//@{
			/** get length */
			inline int length() const {return n_;}
			/** change size */
			inline void resize(int n);
			//@}


			/** {\bf some reasonably efficient math operations} */
			//@{
			/** change sign */
			void negate();

			/** vector addition with result returned through a reference argument */
			inline void add(const DenseSerialVector& other, 
											DenseSerialVector& result) const ;
			/** self-modifying vector addition */
			void add(const DenseSerialVector& other) ;
	
			/** vector subtraction with result returned through a reference argument */
			inline void subtract(const DenseSerialVector& other, 
													 DenseSerialVector& result) const ;
			/** self-modifying vector subtraction */
			void subtract(const DenseSerialVector& other) ;
	
			/** daxpy (z = a*x + y) with result returned through a reference argument */
			inline void daxpy(const DenseSerialVector& other, const TSFReal& a, 
												DenseSerialVector& result) const ;
			/** self-modifying daxpy */
			void daxpy(const DenseSerialVector& other, const TSFReal& a) ;
	
			/** element-by-element multiplication 
			 * with result returned through a reference argument */
			inline void eMult(const DenseSerialVector& other, 
												DenseSerialVector& result) const ;

			/** self-modifying element-by-element multiplication */
			void eMult(const DenseSerialVector& other) ;


            /** absolute value of each element */
            void abs();

            /** return the value of the max element  */
            TSFReal max();

            /** return the value of the min element  */
            TSFReal min();

            /** compute the matlab style ".*" operation, i.e.,
             *      this[i] = y[i] * z[i]  */
            void dotStar(const DenseSerialVector& y, 
                         const DenseSerialVector& z);


            /** compute the matlab style "./" operation, i.e.,
             *      this[i] = y[i] / z[i]  */
            void dotSlash(const DenseSerialVector& y, 
                         const DenseSerialVector& z);

	
			/** multiplication by a scalar
			 * with result returned through a reference argument */
			inline void scalarMult(const TSFReal& scalar, 
														 DenseSerialVector& result) const ;
			/** self-modifying multiplication by a scalar */
			void scalarMult(const TSFReal& scalar);

			/** exponentiation by a scalar
			 * with result returned through a reference argument */
			inline void scalarPow(const TSFReal& scalar, 
														DenseSerialVector& result ) const ;
			/** Self-modifying exponentiation by a scalar */
			void scalarPow(const TSFReal& scalar);

			/** dot product */
			TSFReal dot(const DenseSerialVector& other) const ;

			/** dot product with self */
			TSFReal norm2Squared() const ;

			/** 2-norm */
			TSFReal norm2() const {return sqrt(norm2Squared());}
			
			/** sum elements */
			inline TSFReal sumElements() const ;

			/** return maximum element value */
			TSFReal maxNorm() const ;

			/** set all elements to zero */
			void zero() {setScalar(0.0);}	

			/** set all elements to the given value */
			void setScalar(const TSFReal& a);
			//@}

			/** {\bf overloaded math operators; for performance reasons, 
			 * avoid these in performance-critical
			 * code} */
			//@{
			/** unary minus */
			inline DenseSerialVector operator-() const ;
			/** reflexive addition */
			inline DenseSerialVector& operator+=(const DenseSerialVector& other);
			/** reflexive subtraction */
			inline DenseSerialVector& operator-=(const DenseSerialVector& other);
			/** reflexive scalar mult */
			inline DenseSerialVector& operator*=(const TSFReal& scalar);
			/** reflexive scalar division */
			inline DenseSerialVector& operator/=(const TSFReal& scalar);

			/** addition */
			inline DenseSerialVector operator+(const DenseSerialVector& other) const ;
			/** subtraction */
			inline DenseSerialVector operator-(const DenseSerialVector& other) const ;
			/** dot product */
			inline TSFReal operator*(const DenseSerialVector& other) const ;
			/** scalar mult */
			inline DenseSerialVector operator*(const TSFReal& scalar) const ;
			/** scalar division */
			inline DenseSerialVector operator/(const TSFReal& scalar) const ;
			//@}

			string toString() const ;

			/** write a brief description to string */
			string summary() const ;

			static bool unitTest();

			static TSFTimer blasTimer_;
		private:
			void boundsCheck(int i) const ; 
			TSFReal *x_;
			int n_;

			/* define some utility variables for simplifying calls to fortran BLAS */
			static int one_; 
			static TSFReal onePointZero_;
			static TSFReal negativeOnePointZero_;
		};





	inline DenseSerialVector::DenseSerialVector(int n) : x_(0), n_(n)
		{
			x_ = new TSFReal [n_];
			zero();
		}

	inline DenseSerialVector::DenseSerialVector(int n, TSFReal value) 
		: x_(0), n_(n)
		{
			x_ = new TSFReal [n_];
			for (int i=0; i<n_; i++) x_[i] = value;
		}

	inline DenseSerialVector::DenseSerialVector(int n, const TSFReal* cArray) 
		: x_(0), n_(n)
		{
			x_ = new TSFReal [n_];
			memcpy(x_, (const void*) cArray, (size_t) (n_*sizeof(TSFReal)));
		}

	inline DenseSerialVector::DenseSerialVector(const DenseSerialVector& other) 
		: x_(0), n_(other.n_)
		{
			x_ = new TSFReal [n_];
			memcpy(x_, (const void*) other.x_, (size_t) (n_*sizeof(TSFReal)));
		}

	inline DenseSerialVector::~DenseSerialVector()
		{
			if (x_) delete [] x_;
		}

	inline void DenseSerialVector::resize(int n)
		{
			if (n == n_) return;
			n_ = n;
			if (x_) delete [] x_;
			x_ = new TSFReal[n_];
		}

	inline void DenseSerialVector::boundsCheck(int i) const 
		{
			if (i<0 || i>=n_) TSFError::raise("DenseSerialVector::index out of bounds");
		}

	inline TSFReal& DenseSerialVector::operator[](int i)
		{
#if HAVE_DENSE_VECTOR_BOUNDSCHECK
			boundsCheck(i);
#endif
			return x_[i];
		}

	inline const TSFReal& DenseSerialVector::operator[](int i) const 
		{
#if HAVE_DENSE_VECTOR_BOUNDSCHECK
			boundsCheck(i);
#endif
			return x_[i];
		}

	inline void DenseSerialVector::add(const DenseSerialVector& other, 
																		 DenseSerialVector& result) const 
		{
			result = *this;
			result.add(other);
		}
	
	inline void DenseSerialVector::subtract(const DenseSerialVector& other, 
																					DenseSerialVector& result) const 
		{
			result = *this;
			result.subtract(other);
		}

	inline void DenseSerialVector::daxpy(const DenseSerialVector& other, 
																			 const TSFReal& a, 
																			 DenseSerialVector& result) const 
		{
			result = *this;
			result.daxpy(other, a);
		}

	inline void DenseSerialVector::eMult(const DenseSerialVector& other, 
																			 DenseSerialVector& result) const 
		{
			result = *this;
			result.eMult(other);
		}

	inline void DenseSerialVector::scalarMult(const TSFReal& a, 
																						DenseSerialVector& result) const 
		{
			result = *this;
			result.scalarMult(a);
		}


	inline DenseSerialVector DenseSerialVector::operator-() const 
		{
			DenseSerialVector rtn = *this;
			rtn.negate();
			return rtn;
		}

	inline DenseSerialVector& DenseSerialVector::operator+=(const DenseSerialVector& other) 
	{
		add(other);
		return *this;
	}

	inline DenseSerialVector DenseSerialVector::operator+(const DenseSerialVector& other) const 
		{
			DenseSerialVector rtn = *this;
			rtn.add(other);
			return rtn;
		}

	inline DenseSerialVector& DenseSerialVector::operator-=(const DenseSerialVector& other) 
	{
		subtract(other);
		return *this;
	}

	inline DenseSerialVector DenseSerialVector::operator-(const DenseSerialVector& other) const 
		{
			DenseSerialVector rtn = *this;
			rtn.subtract(other);
			return rtn;
		}

	inline DenseSerialVector& DenseSerialVector::operator*=(const TSFReal& scalar)
	{
		scalarMult(scalar);
		return *this;
	}

	inline DenseSerialVector& DenseSerialVector::operator/=(const TSFReal& scalar)
	{
		if (scalar==0) TSFError::raise("zero divide in DenseSerialVector::operator/=");
		scalarMult(1.0/scalar);
		return *this;
	}

	inline DenseSerialVector DenseSerialVector::operator*(const TSFReal& scalar) const 
		{
			DenseSerialVector rtn = *this;
			rtn.scalarMult(scalar);
			return rtn;
		}

	inline TSFReal DenseSerialVector::operator*(const DenseSerialVector& other) const 
		{
			return dot(other);
		}

	inline DenseSerialVector DenseSerialVector::operator/(const TSFReal& scalar) const 
		{
			DenseSerialVector rtn = *this;
			if (scalar==0) TSFError::raise("DenseSerialVector::operator/ divide by zero");
			rtn.scalarMult(1.0/scalar);
			return rtn;
		}

	inline DenseSerialVector operator*(const TSFReal& scalar, 
																		 const DenseSerialVector& v)
		{
			return v*scalar;
		}

	ostream& operator<<(ostream& os, const DenseSerialVector& v);

}

	


#endif
