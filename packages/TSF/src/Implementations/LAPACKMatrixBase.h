#ifndef LAPACKMATRIXBASE_H
#define LAPACKMATRIXBASE_H

#include "TSFConfig.h"
#include "TSFDenseMatrix.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * Base class for LAPACK matrices.
	 */

	class LAPACKMatrixBase public TSFDenseMatrix
		{
		public:
			/** empty ctor, for when we don't know the matrix size */ 
			LAPACKMatrixBase();

			/** the usual virtual dtor */
			virtual ~LAPACKMatrixBase();
			
			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const = 0 ;

			/** apply inverse operator to a vector in the range space, returning
			 * its preimage as a vector in the domain space. The solve is done
			 * by factoring and backsolving. */
			virtual void applyInverse(const TSFVector& in,
																TSFVector& out) const ;

			/** \name matrix configuration interface */
			//@{
			
			/** configure a dense matrix by setting num rows and num columns */
			virtual void configureDense(int nRows, int nCol) ;
			
			//@}

			/** \name matrix loading interface */
			//@{
			
			/** add to selected elements of a row in the matrix */
			virtual void addToRow(int globalRowIndex,
														int nCols,
														const int* globalColumnIndices,
														const TSFReal* a) = 0 ;

			/** set all elements to zero */
			virtual void zero() ;

			//@}

			/** factor */
			virtual void factor() = 0 ;

			
		protected:
			int nRows_;
			
			int nCols_;

			DenseSerialVector data_;
		};

}

#endif
