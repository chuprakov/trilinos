#ifndef BICGSTABSOLVER_H
#define BICGSTABSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"
#include "TSFPreconditionerFactory.h"
#include "TSFTimeMonitor.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup ConcreteLinearSolvers
	 * Representation-independent BICGSTAB solver. 
	 */
	
	class BICGSTABSolver : public TSFLinearSolverBase
		{
		public:
			/** construct a BICGSTAB solver with 
			 * parameters for max iterations and convergence tolerance */
			BICGSTABSolver(const TSFReal& tol = 1.0e-12, int maxIters = 300);

			/** construct a preconditioned BICGSTAB solver, with 
			 * parameters for max iterations and convergence tolerance */
			BICGSTABSolver(const TSFPreconditionerFactory& pf,
										 const TSFReal& tol = 1.0e-12, int maxIters = 300);

			/** TUVD */
			virtual ~BICGSTABSolver();

			/**
			 * Solve the system with the given RHS, returning the solution 
			 * by reference
			 * argument. The return value is true if the solve succeeded, false
			 * if it failed. 
			 */ 
			virtual bool solve(const TSFLinearOperator& op, 
												 const TSFVector& rhs, 
												 TSFVector& soln) const ;


		private:
			
			virtual bool solveUnpreconditioned(const TSFLinearOperator& op, 
																				 const TSFVector& rhs, 
																				 TSFVector& soln) const ;
			/* residual convergence tolerance */
			TSFReal tol_;
			/* maximum number of iterations */
			int maxIters_;
			/* factory to build a preconditioner */
			TSFPreconditionerFactory preconditionerFactory_;
		};

}


#endif
