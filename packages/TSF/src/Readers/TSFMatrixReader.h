#ifndef TSFMATRIXREADER_H
#define TSFMATRIXREADER_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFMatrixReaderBase.h"

namespace TSF
{
	using std::string;
	using std::istream;

	/** \ingroup MatrixReaders
	 * 
	 */

	class TSFMatrixReader
		{
		public: 
			/** construct with a pointer to a concrete type */
			TSFMatrixReader(TSFMatrixReaderBase* ptr);

			/** read a matrix */
			TSFLinearOperator read(const TSFVectorType& vectorType) const ;

		private:
			TSFSmartPtr<TSFMatrixReaderBase> ptr_;
		};

}


#endif
