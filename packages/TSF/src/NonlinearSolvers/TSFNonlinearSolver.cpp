#include "TSFNonlinearSolver.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"


using namespace TSF;


TSFNonlinearSolver::TSFNonlinearSolver()
	: ptr_(0)
{}

TSFNonlinearSolver::TSFNonlinearSolver(TSFNonlinearSolverBase* ptr)
	: ptr_(ptr)
{}

bool TSFNonlinearSolver::solve(const TSFNonlinearOperator& op, 
														const TSFVector& rhs, 
														TSFVector& soln) const
{
	return ptr_->solve(op - rhs, soln);
}

bool TSFNonlinearSolver::solve(const TSFNonlinearOperator& op, 
															 TSFVector& soln) const
{
	return ptr_->solve(op, soln);
}



bool TSFNonlinearSolver::isNull() const
{
	return ptr_.isNull();
}
