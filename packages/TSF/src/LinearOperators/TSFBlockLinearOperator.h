#ifndef TSFBLOCKLINEAROPERATOR_H
#define TSFBLOCKLINEAROPERATOR_H

#include "TSFConfig.h"
#include "TSFLinearOperatorBase.h"
#include "TSFArray.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 *
	 */

	class TSFBlockLinearOperator : public TSFLinearOperatorBase
		{
		public:
			/** Construct with the specified domain and range */
			TSFBlockLinearOperator(const TSFVectorSpace& domain,
														 const TSFVectorSpace& range);

			/* the usual virtual dtor */
			virtual ~TSFBlockLinearOperator(){;}

			/** get the number of block rows */
			virtual int numBlockRows() const {return nBlockRows_;}

			/** get the number of block columns */
			virtual int numBlockCols() const {return nBlockCols_;}

			/** get the (i,j)-th submatrix */
			virtual void getBlock(int i, int j, TSFLinearOperator& sub) const ;

			/** set the (i,j)-th submatrix */
			virtual void setBlock(int i, int j, 
														const TSFLinearOperator& sub);

			
			/** say whether we are a block operator */
			bool isBlockOperator() const {return true;}
			
			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** apply adjoint operator to a vector in the domain space, returning
			 * a vector in the range space. The default implementation throws an
			 * exception */
			virtual void applyAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/**
			 * Write to a stream 
			 */
			virtual void print(ostream& os) const ;

		protected:
			int nBlockRows_;
			int nBlockCols_;

			TSFArray<TSFArray<TSFLinearOperator> > sub_;
		};
}

#endif
