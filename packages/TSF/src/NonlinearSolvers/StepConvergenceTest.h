#ifndef STEPCONVERGENCETEST_H
#define STEPCONVERGENCETEST_H

#include "TSFVector.h"
#include "TSFConvergenceTest.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup NonlinearSolvers
	 * 
	 */

	class StepConvergenceTest : public TSFConvergenceTestBase
		{
		public:
			/** */
			StepConvergenceTest(const TSFReal& stepTol);

			/** */
			virtual ~StepConvergenceTest(){;}
			
			/** */
			virtual bool testStep(const TSFVector& deltaX) const ;
		private:
			TSFReal stepTol_;
		};
}
#endif
