#ifndef TSFACCESSIBLEVECTOR_H
#define TSFACCESSIBLEVECTOR_H

#include "TSFConfig.h"
#include <typeinfo>
#include "TSFVectorBase.h"

namespace TSF
{
	enum TSFReductionOp {TSFSumReduce, TSFMaxReduce};

	/** \ingroup VectorSubtypes
	 * TSFAccessibleVector is a base class for vector implementations
	 * with efficient access to chunks of data, for example out of core
	 * vectors as well as distributed or serial vectors. An example of
	 * something that would *not* be a TSFAccessibleVector is a proxy
	 * for a vector sitting on a remote machine.  
	 *
	 * The mathematical methods can be implemented efficiently in terms of
	 * BLAS operations on the chunks. 
	 */

	class TSFAccessibleVector : public TSFVectorBase
		{
		public: 
			/** \name Constructor and Destructors */
			//@{
			/** construct with a given space */
			TSFAccessibleVector(const TSFVectorSpace& space);

			/** the usual virtual dtor */
			virtual ~TSFAccessibleVector(){;}
			//@}

			/** \name mathematical methods */
			//@{
			/** x = a*x + y */
			virtual void axpy(const TSFReal& a, const TSFVector& y) ;

			/** multiply by a scalar */
			virtual void scalarMult(const TSFReal& a) ;

			/** element-by-element multiplication (this = this .* y) */
			virtual void dotStar(const TSFVector& x) ;

			/** element-by-element division (this = this ./ y) */
			virtual void dotSlash(const TSFVector& x) ; 

			/** dot product */
			virtual TSFReal dot(const TSFVector& other) const ;

			/** infinity-norm */
			virtual TSFReal normInf() const ;

			/** one-norm */
			virtual TSFReal norm1() const ;

			/** sum all elements */
			virtual TSFReal sumElements() const ;

			/** set all elements to a scalar value */
			virtual void setScalar(const TSFReal& a) ;

			/** find an extreme value, possibly subject to a constraint */
			virtual TSFReal findExtremeValue(MinOrMax type, TSFGeneralizedIndex& location, 
																			 const TSFReal& tol) const ;

			/** absolute value */
			virtual void abs() ;
	
			//@}

			/** \name generating random vectors */
			//@{
			/** Fill a vector with random elements */
			virtual void randomize(const TSFRandomNumberGenerator& r) ;
			//@}	

			/** \name maintenance methods */
			//@{
			/** copy the data from another vector into this vector */
			virtual void acceptCopyOf(const TSFVector& x) ;
			//@}


		protected:
			/** \name Accessing chunks of elements */
			//@{
			/** Set the chunk iterator back to start */
			virtual void rewindChunkIterator() const = 0 ;

			/** Say whether there are more chunks to be processed */
			virtual bool hasMoreChunks() const = 0 ;

			/** Get the next read-only chunk, returning the size via reference argument */
			virtual const TSFReal* getNextChunk(int& chunkSize) const = 0 ;

			/** Get the next read-write chunk, returning the size via reference argument */
			virtual TSFReal* getNextChunk(int& chunkSize) = 0 ;

			/** Get the next read-only chunk of on-processor elements, 
			 * returning the size via reference argument */
			virtual const TSFReal* getNextLocalChunk(int& chunkSize) const = 0 ;

			/** Restore a chunk after use. This is needed for out-of-core vectors. */
			virtual void restoreChunk(const TSFReal* chunk, int chunkSize) const = 0 ;
			//@}

			/** \name reduction methods */
			/** Reduce a calculation across processors */
			virtual TSFReal reduce(const TSFReal& localValue, TSFReductionOp op) const = 0 ;

		};
};

#endif
