#ifndef TSFCONVERGENCETEST_H
#define TSFCONVERGENCETEST_H

#include "TSFConvergenceTestBase.h"
#include "TSFSmartPtr.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup NonlinearSolvers
	 * 
	 */

	class TSFConvergenceTest
		{
		public:
			/** */
			TSFConvergenceTest(TSFConvergenceTestBase* ptr)
				: ptr_(ptr) {}
			
			/** */
			bool testFunctionValue(const TSFVector& f) const 
				{return ptr_->testFunctionValue(f);}

			/** */
			bool testStep(const TSFVector& deltaX) const 
				{return ptr_->testStep(deltaX);}
		private:
			TSFSmartPtr<TSFConvergenceTestBase> ptr_;
		};
}
#endif
