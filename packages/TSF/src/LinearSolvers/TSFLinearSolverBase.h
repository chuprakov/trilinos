#ifndef TSFLINEARSOLVERBASE_H
#define TSFLINEARSOLVERBASE_H

#include "TSFConfig.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup LinearSolverSubtypes
	 */
	
	class TSFLinearSolverBase
		{
		public:
			/** */
			TSFLinearSolverBase(int verbosity=1) : verbosity_(verbosity) {;}
			/** TUVD */
			virtual ~TSFLinearSolverBase();

			/** 
 			 * Solve the system with the given RHS, returning the solution 
			 * by reference
			 * argument. The return value is true if the solve succeeded, false
			 * if it failed. 
			 */ 
			virtual bool solve(const TSFLinearOperator& op, 
												 const TSFVector& rhs, 
												 TSFVector& soln) const = 0 ;	

			/** */
			void setVerbosityLevel(int v) {verbosity_ = v;}

		protected:
			int verbosity_;
		private:
			
		};

}


#endif
