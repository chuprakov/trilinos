#include "TSFLinearSolver.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFPreconditionerFactory.h"
#include "TSFPreconditionerFactoryBase.h"

using namespace TSF;


TSFLinearSolver::TSFLinearSolver()
	: ptr_(0)
{}

TSFLinearSolver::TSFLinearSolver(TSFLinearSolverBase* ptr)
	: ptr_(ptr)
{}

bool TSFLinearSolver::solve(const TSFLinearOperator& op, 
														const TSFVector& rhs, 
														TSFVector& soln) const
{
	return ptr_->solve(op, rhs, soln);
}

bool TSFLinearSolver::isNull() const
{
	return ptr_.isNull();
}

void TSFLinearSolver::setVerbosityLevel(int v) {ptr_->setVerbosityLevel(v);}

