#ifndef LAPACKGENERALMATRIX_H
#define LAPACKGENERALMATRIX_H

#include "TSFConfig.h"
#include "TSFDenseMatrix.h"
#include "DenseSerialVector.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * Linear operator implemented as a LAPACK dense matrix.
	 */

	class LAPACKGeneralMatrix : public TSFMatrixOperator
		{
		public:
			/** Construct with domain and range spaces, which should be 
			 * DenseSerialVectorSpace objects */
			LAPACKGeneralMatrix(const TSFVectorSpace& domain,
													const TSFVectorSpace& range);

            /** Empty Ctor */
            LAPACKGeneralMatrix(){;}
			/** the usual virtual dtor */
			virtual ~LAPACKGeneralMatrix(){;}

            /* added by ptb  */
            inline TSFReal& operator[](int i) {return data_[i];}
            inline const TSFReal& operator[](int i) const {return data_[i];}
            inline TSFReal& operator()(int i, int j) 
              {return data_[i+nRows_*j];}
            inline const TSFReal& operator()(int i, int j) const 
              {return data_[i+nRows_*j];}
			
			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** apply inverse operator to a vector in the range space, returning
			 * its preimage as a vector in the domain space. The solve is done
			 * by factoring and backsolving. */
			virtual void applyInverse(const TSFVector& in,
																TSFVector& out) const ;
			
			/** apply adjoint operator to a vector in the domain space, returning
			 * a vector in the range space. The default implementation throws an
			 * exception */
			virtual void applyAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/** apply inverse adjoint operator */
			virtual void applyInverseAdjoint(const TSFVector& in,
																			 TSFVector& out) const ;

			/** \name matrix loading interface */
			//@{
			
			/** add to selected elements of a row in the matrix */
			virtual void addToRow(int globalRowIndex,
														int nCols,
														const int* globalColumnIndices,
														const TSFReal* a) ;

			/** */
			virtual void setElement(int i, int j, const TSFReal& aij) ;


			/** set all elements to zero */
			virtual void zero() ;

			//@}

			/** factor with a call to LAPACK's getrf() function */
			virtual void factor() ;

			/** write to a stream */
			virtual void print(ostream& os) const ;
			
		protected:
			/** low-level matrix-vector multiply */
			void mvMult(bool transpose, const TSFVector& in, TSFVector& out) const ;

			/** low-level solve */
			void solve(bool transpose, const TSFVector& in, TSFVector& out) const ;

			DenseSerialVector data_;

			TSFArray<int> iPiv_;

			int nRows_;
			
			int nCols_;
		};

}

#endif
