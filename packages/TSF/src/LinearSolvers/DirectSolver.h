#ifndef DIRECTSOLVER_H
#define DIRECTSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"


namespace TSF
{
	using std::ostream;

	/** \ingroup ConcreteLinearSolvers
	 * Representation-independent direct solver. The solve() method invokes
	 * the applyInverse() method of the operator.
	 */ 
	
	class DirectSolver : public TSFLinearSolverBase
		{
		public:
			/** empty ctor is sufficient */
			DirectSolver();

			/** TUVD */
			virtual ~DirectSolver(){;}

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
		};

}


#endif
