#ifndef TSFLINEARSOLVER_H
#define TSFLINEARSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"
#include "TSFSmartPtr.h"

namespace TSF
{
	

	/** \ingroup Solvers
	 * Representation-independent BICGSTAB solver. This is just for testing,
	 * so there is no provision for preconditioning.
	 */
	
	class TSFLinearSolver
		{
		public:
			/** empty ctor */
			TSFLinearSolver();
			/** construct with a ptr to the base class */
			TSFLinearSolver(TSFLinearSolverBase* ptr);

			/** \name Solve methods */
			//@{
			/** 
 			 * Solve the system with the given RHS, returning the solution 
			 * by reference
			 * argument. The return value is true if the solve succeeded, false
			 * if it failed. 
			 */ 
			bool solve(const TSFLinearOperator& op, 
								 const TSFVector& rhs, 
								 TSFVector& soln) const ;	

			//@}
			/** */
			bool isNull() const ;

			/** */
			void setVerbosityLevel(int v);

		private:
			TSFSmartPtr<TSFLinearSolverBase> ptr_;
		};

}


#endif
