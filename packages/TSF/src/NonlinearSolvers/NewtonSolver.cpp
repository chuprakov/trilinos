#include "NewtonSolver.h"
#include "TSFConvergenceTestBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"
#include "TSFOut.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;

NewtonSolver::NewtonSolver(const TSFLinearSolver& linearSolver,
													 int maxIters,
													 const double& stepTol,
													 const double& funcTol)
	: TSFNonlinearSolverBase(maxIters),
		linearSolver_(linearSolver),
		stepTol_(stepTol),
		funcTol_(funcTol)
{;}

NewtonSolver::~NewtonSolver(){}
		
bool NewtonSolver::solve(const TSFNonlinearOperator& F, 
												 TSFVector& soln) const
{
	return solve(*(F.getPtr()), soln);
}

bool NewtonSolver::solve(const TSFNonlinearOperatorBase& F, 
												 TSFVector& soln) const
{
	TSFVector x0 = soln;
	TSFVector f0;
	TSFVector deltaX;
	TSFLinearOperator J;

	int maxBacktracks = 20;

	for (int i=0; i<maxIters_; i++)
		{
			/* compute the function value */
			F.apply(x0, f0);

			double initialFNorm = f0.norm2();
			if (initialFNorm < funcTol_)
				{
					soln = x0;
					TSFOut::printf("Newton's method converged with func norm %g after %d steps\n",
												 initialFNorm, i++);
					return true;
				}

			/* solve for the full newton step */
			J = F.derivative(x0);
			J.applyInverse(linearSolver_, f0, deltaX);

			double dx = deltaX.norm2();
			TSFOut::printf("iter=%d newton step norm=%g\n",
										 i, dx);
			/* test for step convergence */
			if (dx < stepTol_) 
				{
					soln = x0 - deltaX;
					TSFOut::printf("Newton's method converged with stepsize %g after %d steps\n",
												 dx, i++);
					return true;
				}

			/* The Newton direction is a descent direction of norm(F), but
			 * it's possible that the Newton step has overshot. If the full 
			 * step does not decrease the residual norm, we halve the
			 * step length until a decrease is found. */
			double stepFrac = 1.0;
			bool accepted = false;
			for (int b=0; b<maxBacktracks; b++)
				{
					TSFVector xTry = x0 - stepFrac * deltaX;


					F.apply(xTry, f0);
					double newFNorm = f0.norm2();
					TSFOut::printf("backtrack %d fNew=%g f=%g\n", b, newFNorm,
												 initialFNorm);
					if (newFNorm < initialFNorm) 
						{
							if (newFNorm < funcTol_)
								{
									soln = xTry;
									TSFOut::printf("Newton's method converged with func norm %g after %d steps\n",
												 newFNorm, i++);
									return true;
								}
							x0 = xTry;
							accepted = true;
							break;
						}
					stepFrac *= 0.5;
				}
			if (!accepted) 
				{
					return false;
				}

		}
	
	soln = x0;
	return false;
}
