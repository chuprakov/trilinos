#ifndef TSFLINEARPROBLEM_H
#define TSFLINEARPROBLEM_H

#include "TSFConfig.h"
#include "TSFLinearProblemBase.h"
#include "TSFSmartPtr.h"


namespace TSF
{
	

	/** \ingroup LinearProblem
	 *
	 * TSFLinearProblem is a user-level representation of a linear problem
	 * A x = b, where A is a TSFLinearOperator and x and b are TSFVectors. 
	 * 
	 */

	class TSFLinearProblem
		{
		public:
			/** construct with a pointer to a linear problem subtype */
			TSFLinearProblem(TSFLinearProblemBase* ptr);

			/** returns the right-hand side vector b */
			TSFVector getRHS() const ;

			/** For testing, returns the known solution x to the problem A x = b */
			TSFVector getKnownSolution() const ;

			/** returns the linear operator A */
			TSFLinearOperator getOperator() const ;

		private:
			TSFSmartPtr<TSFLinearProblemBase> ptr_;
		};
}

	
#endif
