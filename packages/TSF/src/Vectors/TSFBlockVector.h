#ifndef TSFBLOCKVECTOR_H
#define TSFBLOCKVECTOR_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

namespace TSF
{
	/** \ingroup VectorSubtypes
	 * Block vectors
	 */

	class TSFBlockVector : public TSFVectorBase
		{
		public: 
			/** construct with a given space */
			TSFBlockVector(const TSFVectorSpace& space);

			/** the usual virtual dtor */
			virtual ~TSFBlockVector(){;}
	
			/** \name Block access */
			//@{
			/** return the number of subvector blocks */
			virtual int numBlocks() const ;

			/** return the i-th subvector */
			virtual void getBlock(int i, const TSFVector& self, 
														TSFVector& sub) const ;

			/** set the i-th subvector */
			virtual void setBlock(int i, const TSFVector& sub);
			//@}

			/** \name mathematical methods */
			//@{
			/** x = a*x + y */
			virtual void axpy(const TSFReal& a, const TSFVector& y) ;

			/** multiply by a scalar */
			virtual void scalarMult(const TSFReal& a) ;

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
			//@}	


		
#if HAVE_RTOP
			/** \name General reduction and transformation operators */
			//@{
			/** apply a transformation operator */
			virtual void apply(RTOp op) ;
			//@}
#endif
					
			/** \name access to elements */
			//@{
			/** read-only access to a single element */
			virtual const TSFReal& getElement(int globalIndex) const ;

			/** read-write access to a single element */
			virtual TSFReal& setElement(int globalIndex) ;

			/** set a block of elements */
			virtual void setElements(int n, const int* globalIndices, 
															 const TSFReal* values) ;

			/** get a block of elements */
			virtual void getElements(int n, const int* globalIndices, 
															 TSFReal* values) const ;
			
			/** add to a block of elements */
			virtual void addToElements(int n, const int* globalIndices, 
																 const TSFReal* values) ;
			//@}	

		

			/** \name generating random vectors */
			//@{
			/** Fill a vector with random elements */
			virtual void randomize(const TSFRandomNumberGenerator& r) ;
			//@}	

			/** \name Hooks for parallel support */
			//@{
			/** gather valid ghost values from other procs. Default is a no-op */
			virtual void synchronizeGhostValues() const ;

			/** mark ghost values as invalid, meaning that they need to be
			 * synchronized. Default is a no-op.  */
			virtual void invalidateGhostValues() ;
			//@}

			/** \name maintenance methods */
			//@{
			/** virtual copy ctor */
			virtual TSFVectorBase* deepCopy() const ;

			/** copy the data from another vector into this vector */
			virtual void acceptCopyOf(const TSFVector& x) ;

			/** write to stream */
			virtual ostream& print(ostream& os) const ;
			//@}
		protected:
			/* array of subvectors */
			TSFArray<TSFVector> subvectors_;
		};
};

#endif
