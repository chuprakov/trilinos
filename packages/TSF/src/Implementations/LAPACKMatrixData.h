#ifndef LAPACKMATRIXDATA_H
#define LAPACKMATRIXDATA_H

#include "TSFConfig.h"
#include "DenseSerialVector.h"
#include "TSFSmartPtr.h"
#include "TSFArray.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * Base class for general LAPACK matrices.
	 */

	class LAPACKGeneralMatrixData
		{
		public:
			/** empty ctor */
			LAPACKGeneralMatrixData();

			/** make a deep copy of the data */
			LAPACKGeneralMatrixData deepCopy() const ;

			/** number of rows */
			int nRows_;
			
			/** number of columns */
			int nCols_;

			/** the data */
			TSFSmartPtr<DenseSerialVector> data_;

			/** pivot vector */
			TSFSmartPtr<TSFArray<int> > iPiv_;
		};

}

#endif
