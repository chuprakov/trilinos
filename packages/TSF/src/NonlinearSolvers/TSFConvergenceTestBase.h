#ifndef TSFCONVERGENCETESTBASE_H
#define TSFCONVERGENCETESTBASE_H

#include "TSFConfig.h"
#include "TSFVector.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup NonlinearSolvers
	 * 
	 */

	class TSFConvergenceTestBase
		{
		public:
			/** */
			TSFConvergenceTestBase(){;}

			/** */
			virtual ~TSFConvergenceTestBase(){;}
			
			/** */
			virtual bool testFunctionValue(const TSFVector& /* f */) const 
				{return false;}

			/** */
			virtual bool testStep(const TSFVector& /* deltaX */) 
				const {return false;}
		private:
		};
}
#endif
