#ifndef TSFMATRIXREADERBASE_H
#define TSFMATRIXREADERBASE_H

#include "TSFConfig.h"
#include <string>
#include <iostream>
#include "TSFLinearOperator.h"
#include "TSFVectorType.h"
#include "TSFVectorSpaceBase.h"

namespace TSF
{
	using std::string;
	using std::istream;

	/** \ingroup MatrixReaders
	 * */

	class TSFMatrixReaderBase
		{
		public: 
			/** construct with the name of the file that stores the matrix */
			TSFMatrixReaderBase(const string& filename) : filename_(filename) {;}

			/** virtual dtor */
			virtual ~TSFMatrixReaderBase(){;}

			/** read a matrix */
			virtual TSFLinearOperator read(const TSFVectorType& vectorType) const = 0 ;

		protected:
			string filename_;
		};
}


#endif
