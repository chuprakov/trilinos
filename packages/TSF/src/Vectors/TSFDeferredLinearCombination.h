#ifndef TSFDEFERREDLINEARCOMBINATION_H
#define TSFDEFERREDLINEARCOMBINATION_H

#include "TSFVector.h"
#include "TSFSmartPtr.h"
#include "TSFTimeMonitor.h"

namespace TSF
{
		
	

	/** \ingroup Vector 
	 * Stores a linear combination of vectors for subsequent evaluation. This
	 * is used to avoid copies and temporary vectors during the evaluation of
	 * overloaded vector operations.
	 */
	
	class TSFDeferredLinearCombination
		{
		public:
			TSFDeferredLinearCombination(const TSFReal& coeff, const TSFVector& v);

			/** unary minus, changes signs of all terms  */
			TSFDeferredLinearCombination operator-() const ;

			/** adds two sums */
			TSFDeferredLinearCombination operator+(const TSFDeferredLinearCombination& op) const ;
			
			/** subtracts two sums */
			TSFDeferredLinearCombination operator-(const TSFDeferredLinearCombination& op) const ;

			/** adds a new vector into the sum */
			TSFDeferredLinearCombination operator+(const TSFVector& v) const ;
			
			/** subtracts a vector from the sum */
			TSFDeferredLinearCombination operator-(const TSFVector& v) const ;

			/** multiplies the sum by a scalar constant */
			TSFDeferredLinearCombination operator*(const TSFReal& a) const ;

			/** takes a dot product with another sum */
			TSFReal operator*(const TSFDeferredLinearCombination& b) const ;

			/** takes a dot product with a vector */
			TSFReal operator*(const TSFVector& b) const ;

			/** add a new term to the list */
			void add(const TSFReal& coeff, const TSFVector& v);

			/** carry out the operation, returning the result by ref argument */
			void evaluateIntoExistingVector(TSFVector& result) const ;

			/** carry out the operation, returning the result as a new vector */
			TSFVector evaluateIntoNewVector() const ;

			/**  2-norm squared (vector dotted with self) */
			TSFReal norm2Squared() const ;

			/** 2-norm */
			TSFReal norm2() const ;

			/** 1-norm */
			TSFReal norm1() const ;

			/** inf-norm */
			TSFReal normInf() const ;
			
		private:
			TSFArray<bool> markLHSVectors(const TSFVector& x) const ;
			
			TSFArray<TSFReal> coeffs_;
			TSFArray<TSFVector> vectors_;
		};

	/** \relates TSFVector */
	inline TSFDeferredLinearCombination 
		operator+(const TSFVector& v, 
							const TSFDeferredLinearCombination& op)
		{return op+v;}

	/** \relates TSFVector */
	inline TSFDeferredLinearCombination 
		operator-(const TSFVector& v, 
							const TSFDeferredLinearCombination& op)
		{return -op + v;}

	/** \relates TSFVector */
	inline TSFDeferredLinearCombination 
		operator*(const TSFReal& a, 
							const TSFDeferredLinearCombination& op)
		{return op*a;}

	/** \relates TSFVector */
	TSFDeferredLinearCombination operator+(const TSFVector& a, 
							const TSFVector& b);
	
	/** \relates TSFVector */
	TSFDeferredLinearCombination operator-(const TSFVector& a, 
																 const TSFVector& b);

	/** \relates TSFVector */
	inline TSFReal operator*(const TSFVector& a, 
													const TSFDeferredLinearCombination& b)
		{return b*a;}

	/** \relates TSFVector */
	inline TSFReal operator*(const TSFVector& a, const TSFVector& b)
		{return a.dot(b);}

	/** \relates TSFVector */
	inline TSFDeferredLinearCombination operator*(const double& a,
																								const TSFVector& x)
		{return TSFDeferredLinearCombination(a, x);}
	
	/** \relates TSFVector */
	inline TSFDeferredLinearCombination operator*(const TSFVector& x,
																								const double& a)
		{return a*x;}
}

#endif
