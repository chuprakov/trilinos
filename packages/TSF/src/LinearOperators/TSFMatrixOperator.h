#ifndef TSFMATRIXOPERATOR_H
#define TSFMATRIXOPERATOR_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFNonDupArray.h"
#include "TSFLinearOperatorBase.h"

namespace TSF
{
	class TSFPreconditioner;
	

	/** \ingroup LinearOperatorSubtypes
	 * Base class for operators that are implemented in terms of a matrix.
	 */

	class TSFMatrixOperator : public TSFLinearOperatorBase
		{
		public:
			/** construct with domain and range spaces */
			TSFMatrixOperator(const TSFVectorSpace& domain, 
												const TSFVectorSpace& range);

			/** the usual virtual dtor */
			virtual ~TSFMatrixOperator();

			/** identify self as a matrix operator */
			virtual bool isMatrixOperator() const {return true;}

			/** */
			const TSFSmartPtr<const TSFMatrixOperator> getMatrix() const {return TSFSmartPtr<const TSFMatrixOperator>(this,false);}

			/** \name matrix configuration interface */
			//@{
			/** inform caller if a full graph is needed to configure this matrix */
			virtual bool requiresGraph() const {return false;}

			/** inform caller if a bandwidth array
					is needed to configure this matrix */
			virtual bool requiresBandwidth() const {return false;}

			/** inform caller if this matrix type can handle non-square matrices */
			virtual bool supportsNonSquare() const {return false;}

			/** Set the full sparsity graph, giving the column indices for each
			 * row owned by the current processor. */
			virtual void setGraph(int nLocalRows, const int* bandwidth,
														const int** columnIndices) ;

			/** set the columns to be used in a given row */
			virtual void setRowStructure(int globalRowIndex, int bandwidth,
																	 const int* columnIndices);

			/** set the bandwith of each row */
			virtual void setBandwidth(int nLocalRows, const int* bandwidth) ;

			/** hook for implementation-dependent structure finalization call */
			virtual void freezeStructure() ;

			/** hook for implementation-dependent structure finalization call */
			virtual void freezeValues() ;
			//@}

			/** \name matrix loading interface */
			//@{
			
			/** add to selected elements of a row in the matrix */
			virtual void addToRow(int globalRowIndex,
														int nCols,
														const int* globalColumnIndices,
														const TSFReal* a) = 0 ;

			/** set a single element */
			virtual void setElement(int i, int j, const TSFReal& aij) ;

			/** set all elements to zero */
			virtual void zero() = 0 ;
			//@}

			/** \name incomplete factorization preconditioning interface */
			//@{
			/** create a k-level incomplete factorization. Default is to throw
			 * an error. */
			virtual void getILUKPreconditioner(int fillLevels,
																				 int overlapFill,
																				 TSFPreconditioner& rtn) const ;
			//@}

			/** \name factoring interface */
			//@{
			/** indicate whether the matrix is now stored in factored form */
			virtual bool isFactored() const ;
			
			/** factor the matrix */
			virtual void factor() ;
			//@}
			
		protected:
			bool isFactored_;

		private:
			TSFMatrixOperator(); // Not defined and not to be called!
			
		};

}

#endif
