#ifndef TSFTESTPROBLEMBASE_H
#define TSFTESTPROBLEMBASE_H

#include "TSFConfig.h"
#include "TSFLinearProblem.h"

namespace TSF
{
	

	/** \ingroup TestProblems
	 * TSFTestProblem is the base class for linear problems. For a test
	 * problem, one must supply a known solution vector in addition to
	 * the operator and RHS. 
	 */

	class TSFTestProblemBase
		{
		public:
			/** empty ctor only */
			TSFTestProblemBase(){;}

			/** virtual dtor */
			virtual ~TSFTestProblemBase(){;}

			/** return the known solution of the test problem */
			virtual TSFVector getKnownSolution() = 0 ;

		private:
		};
}

