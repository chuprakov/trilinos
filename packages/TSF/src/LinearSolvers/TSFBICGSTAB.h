#ifndef BICGSTABSOLVER_H
#define BICGSTABSOLVER_H

#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

namespace TSF
{
	/** \ingroup Hilbert
	 * Representation-independent BICGSTAB solver. This is just for testing,
	 * so there is no provision for preconditioning.
	 */
	
	class BICGSTAB
		{
		public:
			/** construct a BICGSTAB solver to work with a given operator, with 
			 * parameters for max iterations and convergence tolerance */
			BICGSTAB(const TSFLinearOperator& op,
							 const double& tol = 1.0e-12, int maxIters = 300);
			
			/**
			 * Solve the system with the given RHS, returning the solution 
			 * by reference
			 * argument. The return value is true if the solve succeeded, false
			 * if it failed. 
			 */ 
			bool solve(const TSFVector& rhs, 
								 TSFVector& soln) const ;
		private:
			TSFLinearOperator op_;
			double tol_;
			int maxIters_;
		};

}


#endif
