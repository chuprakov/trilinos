#ifndef FUNCTIONVALUECONVERGENCETEST_H
#define FUNCTIONVALUECONVERGENCETEST_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFConvergenceTestBase.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup NonlinearSolvers
	 * 
	 */

	class FunctionValueConvergenceTest : public TSFConvergenceTestBase
		{
		public:
			/** */
			FunctionValueConvergenceTest(const TSFReal& fTol);

			/** */
			virtual ~FunctionValueConvergenceTest(){;}
			
			/** */
			virtual bool testFunctionValue(const TSFVector& f) const ;
		private:
			TSFReal fTol_;
		};
}
#endif
