#include "TSFNonlinearSolverBase.h"
#include "TSFError.h"

using namespace TSF;

TSFNonlinearSolverBase::TSFNonlinearSolverBase(int maxIters)
	: maxIters_(maxIters)
{
	if (maxIters_ <= 0) TSFError::raise("TSFNonlinearSolverBase given non-positive max iteration count");
}
