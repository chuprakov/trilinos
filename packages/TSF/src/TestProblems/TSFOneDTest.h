#ifndef TSFONEDTESTPROBLEM_H
#define TSFONEDTESTPROBLEM_H

#include "TSFConfig.h"
#include "TSFDefaultMatrixProblem.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup CoreSubtypes
	 * TSFMatrixProblem
	 * 
	 */

	class TSFOneDTestProblem : public TSFDefaultMatrixProblem
		{
		public:
			/** empty ctor only */
			TSFOneDTestProblem(TSFReal a, TSFReal b, int n, 
												 TSFMatrixOperator* matrix);

			/** virtual dtor */
			virtual ~TSFOneDTestProblem(){;}

		protected:
			/** */
			virtual int nGlobalRows() const {return nGlobal_;}

			/** */
			virtual int nLocalRows() const {return nLocal_;}

			/** */
			virtual int lowestLocalRow() const {return myLowestRow_;}

			TSFReal h_;

			TSFReal a_;

			TSFReal b_;
			
			int nGlobal_;

			int nLocal_;

			int myLowestRow_;
		};
}

#endif
