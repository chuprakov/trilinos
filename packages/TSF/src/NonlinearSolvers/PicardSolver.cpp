#include "PicardSolver.h"
#include "TSFConvergenceTestBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;

PicardSolver::PicardSolver(int maxIters,
													 const double& tol,
													 const double& relax)
	: TSFNonlinearSolverBase(maxIters),
		tol_(tol),
		relax_(relax)
{;}

PicardSolver::~PicardSolver(){;}		

bool PicardSolver::solve(const TSFNonlinearOperator& F, 
												 TSFVector& soln) const
{

	TSFVector fVal;

	for (int i=0; i<maxIters_; i++)
		{
			TSFVector x0 = soln;
			F.apply(x0, fVal);
			soln = relax_ * x0 + (1.0-relax_)*fVal;
			double res = (soln-x0).norm2();
			if (res < tol_) return true;
		}

	return false;
}

bool PicardSolver::solve(const TSFNonlinearOperatorBase& F, 
												 TSFVector& soln) const
{

	TSFVector fVal;

	for (int i=0; i<maxIters_; i++)
		{
			TSFVector x0 = soln;
			F.apply(x0, fVal);
			soln = relax_ * x0 + (1.0-relax_)*fVal;
			double res = (soln-x0).norm2();
			if (res < tol_) return true;
		}

	return false;
}
