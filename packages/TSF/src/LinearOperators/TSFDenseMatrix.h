#ifndef TSFDENSEMATRIX_H
#define TSFDENSEMATRIX_H

#include "TSFConfig.h"
#include "TSFMatrixOperator.h"
#include "TSFLinearOperator.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * Base class for operators that are implemented in terms of a dense matrix.
	 */

	class TSFDenseMatrix : public TSFMatrixOperator
		{
		public:
			/** empty ctor, for when we don't know the matrix size */ 
			TSFDenseMatrix();

			/** the usual virtual dtor */
			virtual ~TSFDenseMatrix(){;}

			/** \name matrix configuration interface */
			//@{
			
			/** inform caller if the matrix is dense */
			virtual bool isDense() const {return true;}

			/** configure a dense matrix by setting num rows and num columns */
			virtual void configureDense(int nRows, int nCol) = 0 ;
			
			/** set an element */
			virtual void setElement(int i, int j, const TSFReal& aij) = 0 ;

			//@}

			/** \name factoring interface */
			//@{
			/** */
			virtual bool isFactored() const {return isFactored_;}


			
			/** */
			virtual void factor() = 0 ;
			//@}
			
		protected:
			/** flag indicating whether the matrix has been factored */
			bool isFactored_;
		};

}

#endif
