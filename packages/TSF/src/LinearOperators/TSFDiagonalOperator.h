#ifndef TSFDIAGONALOPERATOR_H
#define TSFDIAGONALOPERATOR_H

#include "TSFConfig.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * A TSFDiagonalOperator 
	 */

	class TSFDiagonalOperator : public TSFLinearOperatorBase
		{
		public:
			/** */
			TSFDiagonalOperator(const TSFVector& diagonalValues);
										 
			/** the usual virtual dtor */
			virtual ~TSFDiagonalOperator(){;}

			/** apply does an element-by-element multiply between the input vector
			 * and the diagonal values. */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** applyAdjoint does an element-by-element multiply between 
			 * the input vector and the diagonal values. */
			virtual void applyAdjoint(const TSFVector& in, 
																TSFVector& out) const ;

			/** applyInverse does an element-by-element multiply between 
			 * the input vector and the reciprocals of the diagonal values. */
			virtual void applyInverse(const TSFVector& in, 
																TSFVector& out) const ;
			/**
			 * Write to a stream 
			 */
			virtual void print(ostream& os) const ;
		protected:
			TSFVector diagonalValues_;
		};
}

#endif
