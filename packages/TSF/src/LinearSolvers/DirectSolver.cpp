#include "DirectSolver.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFBlas.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"

using namespace TSF;


DirectSolver::DirectSolver()
	: TSFLinearSolverBase()
{;}

bool DirectSolver::solve(const TSFLinearOperator& op,
												 const TSFVector& b,
												 TSFVector& soln) const 
{
	op.applyInverse(b, soln);
	return true;
}


