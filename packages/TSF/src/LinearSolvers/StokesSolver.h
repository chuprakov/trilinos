#ifndef STOKESSOLVER_H
#define STOKESSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"
#include "TSFLinearSolver.h"
#include "TSFPreconditionerFactory.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup ConcreteLinearSolvers
	 */

	class StokesSolver : public TSFLinearSolverBase
		{
		public:
			StokesSolver(const TSFLinearSolver& innerSolver,
														const TSFLinearSolver& outerSolver);

			virtual ~StokesSolver(){;}

			/**
			 * Solve the system with the given RHS, returning the solution 
			 * by reference
			 * argument. The return value is true if the solve succeeded, false
			 * if it failed. The input operator is expected to be a 2x2 block operator,
			 * and the input rhs is expected to be a
			 */ 
			virtual bool solve(const TSFLinearOperator& op, 
												 const TSFVector& rhs, 
												 TSFVector& soln) const ;
		private:
			
			TSFLinearSolver innerSolver_;

			TSFLinearSolver outerSolver_;
		};

}


#endif
