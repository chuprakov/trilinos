#ifndef TSFVECTOR_H
#define TSFVECTOR_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFTimer.h"
#include "TSFTimeMonitor.h"
#include "TSFRandomNumberGenerator.h"
#include "SystemRand.h"
#include "TSFDeferredCopy.h"
#include "TSFArray.h"
#include "TSFVectorBase.h"
#include <iostream>
#include <string>

namespace TSF
{
	class DenseSerialVector;
	class TSFDeferredLinearCombination;
	
	using std::string;
	using std::ostream;
	

	/** \ingroup Vector
	 * TSFVector:  User-level handle class for vectors. 
	 *
	 * <b> Constructing a TSFVector </b>
	 * 
	 * Ordinarily, you will never construct a TSFVector directly from a derived type. 
	 * Rather, the createMember() method of TSFVectorSpace is used to build a vector
	 * of the appropriate type, for example,
	 * \code
	 * TSFVector x = mySpace.createMember();
	 * \endcode
	 * This hides from you all the ugly details of creating a particular concrete type.
	 *
	 * You will frequently create an empty vector to be filled in later, for example,
	 * \code 
	 * TSFVector y;
	 * \endcode
	 * Note that this vector isn't just empty, it's null. Not only does it have no
	 * values assigned, it does not have a concrete type. An attempt to set an element
	 * of or operate on a null vector will result in an error. What you can do is
	 * assign another vector to it,
	 * \code
	 * TSFVector y;
	 * y = x.deepCopy();
	 * \endcode
	 * or fill it by passing it as a return-by-reference argument to a mathematical
	 * operation,
	 * \code
	 * TSFVector y;
	 * x.scalarMult(3.14159, y);
	 * \endcode
	 * 
	 *
	 * Another way a TSFVector can be created is as the result of an operation between
	 * existing vectors, for example
	 * \code
	 * TSFVector z = a*x + b*y;
	 * \endcode
	 * See the section on operator overloading for more information on this. 
	 * 
	 *
	 * <b> Element access </b>
	 *
	 * There are a number of methods to get or set the value of an element or group 
	 * of elements at a particular index or group of indices. In a distributed environment,
	 * there is a distinction between an element's <i> local </i> index and its
	 * <i> global </i> index. All the element access methods of TSFVector use global
	 * indices. 
	 *
	 * Element access methods must make virtual function calls, so they
	 * should be used as sparingly as possible. In particular, you
	 * should <i> never </i> do mathematical operations on the whole
	 * vector using the element access methods. Use one of the
	 * predefined operations (e.g., eMult, update) or a custom reduction
	 * and transformation operator instead. If you are needing to set
	 * elements, such as when creating the load vector in a
	 * finite-element problem, it is more efficient to set groups of
	 * elements rather than single elements.
	 *
	 * <b> Overloaded operators </b>
	 *
	 * TSFVector supports the full set of overloaded operators for vector and scalar-vector
	 * arithmetic. Overloaded operators are not often used in scientific computing
	 * because a naive implementation incurs a large performance penalty due to 
	 * the creation and copying of temporary vectors. TSFVector operators use a design
	 * which avoids the creation of such temporary vectors by flattening the expression
	 * tree. Small temporary objects are created in the process, so for small problems
	 * (under \f$N \approx 1000 \f$) there is still a noticeable performance penatly. 
	 * However, for larger problems the performance penalty is less than a few percent
	 * compared to FORTRAN BLAS, and decreases quickly with problem size. 
	 * */

	class TSFVector
		{
		public: 
			/** \name Constructors, Destructors, and Assignment Operators */
			//@{
			/** empty ctor. Constructs a null vector */
			TSFVector();

			/** Construct with a pointer to a derived vector type. */
			TSFVector(TSFVectorBase* ptr);
			
			/** construct with the result of an operation */
			TSFVector(const TSFDeferredLinearCombination& op);

			/** construct with the result of an operation */
			TSFVector(const TSFDeferredCopy& copy);
			
			/** assign the result of an operation into this vector */
			TSFVector& operator=(const TSFDeferredLinearCombination& op);

			/** assign the result of an operation into this vector */
			TSFVector& operator=(const TSFDeferredCopy& copy);
			//@}


			/** \name Vector space support */
			//@{
			/** return the vector space in which this vector lives. */
			const TSFVectorSpace& space() const ;

			/** indicate whether this vector is a member of a given space */
			bool isMemberOf(const TSFVectorSpace& space) const ;
			//@}

			/** \name Block access */
			//@{
			/** return the number of subvector blocks */
			int numBlocks() const ;

			/** return the i-th subvector */
			TSFVector getBlock(int i) const ;

			/** set the i-th subvector */
			TSFVector& setBlock(int i, const TSFVector& sub);
			
			//@}


			/** \name Setting, getting, and adding to individual elements or groups of elements. */
			//@{
			/**  Read-only access to an element specified by its global index  */
			const TSFReal& operator[](int globalIndex) const ;

			/** Read-write access to an element specified by its global index  */
			TSFReal& operator[](int globalIndex) ;
			
			/** pack selected elements into a vector */
			void getElements(const TSFArray<int>& globalIndices, 
											 DenseSerialVector& sub) const ;
			
			/** set several elements */
			void setElements(const TSFArray<int>& globalIndices, 
											 const DenseSerialVector& sub);
			
			/** add to several elements */
			void addToElements(const TSFArray<int>& globalIndices, 
												 const DenseSerialVector& sub);

			/** pack selected elements into a vector */
			void getElements(const int* globalIndices, const int length,
											 DenseSerialVector& sub) const ;
			
			/** set several elements */
			void setElements(const int* globalIndices,
											 const DenseSerialVector& sub);
			
			/** add to several elements */
			void addToElements(const int* globalIndices,
												 const DenseSerialVector& sub);

			/** add to a selected element */
			void addToElement(int globalIndex, const TSFReal& val);
			//@}
			
			/** \name Math operations */
			//@{
			/** axpy (this = a*x + y) */
			inline void axpy(const TSFReal& a, const TSFVector& x, const TSFVector& y)
				{acceptCopyOf(y); selfModifyingAxpy(a, x);}

			/** multiplication by a scalar (this = a*x) */
			inline void scalarMult(const TSFReal& a, const TSFVector& x) 
				{acceptCopyOf(x); selfModifyingScalarMult(a);}
			
			/** addition (this = x + y) */
			inline void add(const TSFVector& x, const TSFVector& y) {axpy(1.0, x, y);}

			/** subtraction (this = x - y) */
			inline void subtract(const TSFVector& x, const TSFVector& y) {axpy(-1.0, y, x);}

			/** dot product with another vector */
			TSFReal dot(const TSFVector& other) const ;
	
			/** 2-norm */
			TSFReal norm2() const ;

			/** 1-norm */
			TSFReal norm1() const ;

			/** inf-norm */
			TSFReal normInf() const ;
			
			/** set all elements to zero */
			void zero();	

			/** sum all elements */
			TSFReal sumElements() const ;

			/** set all elements to a scalar value */
			void setScalar(const TSFReal& a);

			//@}

			/** \name generating random vectors */
			//@{
			/** Fill a vector with random elements */
			void randomize(const TSFRandomNumberGenerator& r = new SystemRand()) ;
			//@}
			
			/** \name copying */
			//@{
			/** make a deep copy of this vector. The execution of the copy is
			 * deferred until assignment, so that we can test whether we need
			 * to allocate space for the copy. */
			TSFDeferredCopy copy() const ;
			//@}

			/** \name output */
			//@{
			/** print, called by stream output */
			void print(ostream& os) const ;

			/** write a representation of this object to a string */
			string toString() const ;
			//@}


			/** \name introspection */
			//@{
			/** determine if a vector has not been initialized */
			bool isNull() const ;

			/** see if two vectors share the same pointer */
			bool isIdenticalTo(const TSFVector& other) const ;
			//@}

			/** \name Hooks for parallel support */
			//@{
			/** gather valid ghost values from other procs */
			void synchronizeGhostValues() const ;

			/** mark ghost values as invalid, meaning that they need to be
			 * synchronized */
			void invalidateGhostValues() ;
			//@}

			/** \name dangerous developer-only methods. */
			//@{
			/** read-only pointer access */
			const TSFSmartPtr<TSFVectorBase>& smartPtr() const {return ptr_;}
			/** read-write pointer access */
			TSFSmartPtr<TSFVectorBase>& smartPtr() {return ptr_;}
			/** read-only pointer access */
			const TSFVectorBase* ptr() const {return &(*ptr_);}
			/** read-write pointer access */
			TSFVectorBase* ptr() {return &(*ptr_);}
			//@}

			/** \name internal math operations */
			//@{
			/** copy another's contents into self */
			void acceptCopyOf(const TSFVector& x);
			/** set self = self + a*x */
			void selfModifyingAxpy(const TSFReal& a, const TSFVector& x);
			/** set self = a*self */
			void selfModifyingScalarMult(const TSFReal& a);
			//@}
		private:

			TSFSmartPtr<TSFVectorBase> ptr_;

			/** timer for math operations */
			static TSFTimer& opTimer();
			/** timer for deep copies */
			static TSFTimer& copyTimer();
		};

	/** \relates TSFVector
	 * write an TSFVector to an output stream
	 */
	inline ostream& operator<<(ostream& os, const TSFVector& x)
		{
			x.print(os);
			return os;
		}

}


#endif
