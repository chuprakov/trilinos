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
			LAPACKGeneralMatrixData() : is_factored_(false) {}

			/** make a deep copy of the data */
			LAPACKGeneralMatrixData deepCopy() const { assert(0); return *this; }

			/** number of rows */
			int nRows_;
			
			/** number of columns */
			int nCols_;

			/** the data */
			TSFSmartPtr<DenseSerialVector> data_;

			/** the factor's data */
			TSFSmartPtr<DenseSerialVector> factor_data_;

			/** pivot vector */
			TSFSmartPtr<TSFArray<int> > iPiv_;
			
			/** Is factored! */
			bool is_factored_;
			
		};

}

#endif
