#include "TSFMatrixProblem.h"

#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"

using namespace TSF;



TSFMatrixProblem::TSFMatrixProblem(TSFMatrixOperator* matrix)
	: TSFLinearProblemBase(), 
	op_(matrix), 
	matrix_(matrix), 
	localRowIndices_(0),
	graph_(0),
	matrixReady_(false)
{}


void TSFMatrixProblem::buildMatrixStructure() const 
{
	/* if the graph is needed for this matrix type, create it and set it */
	if (matrix_->requiresGraph())
		{
			//
		}
	else if (matrix_->requiresBandwidth())
		{
			TSFSmartPtr<TSFArray<int> > bandwidth = formBandwidth();
			matrix_->setBandwidth((int) bandwidth->size(), &((*bandwidth)[0]));
		}
	
	/* tell the matrix that we're done with configuration */
	matrix_->freezeStructure();

	/* set all elements to zero */
	matrix_->zero();
}


TSFLinearOperator TSFMatrixProblem::getOperator() const 
{
	if(!matrixReady_)
		{
			buildMatrixStructure();
			fillMatrix();
			matrix_->freezeValues();
			matrixReady_ = true;
		}
	return op_;
}

TSFVector TSFMatrixProblem::getRHS() const 
{
	return getRHS(getOperator().range());
}


TSFVector TSFMatrixProblem::getRHS(const TSFVectorSpace& space) const 
{
	TSFVector rhs = space.createMember();
	fillRHS(rhs);
	return rhs;
}


TSFVector TSFMatrixProblem::getKnownSolution() const 
{
	return getKnownSolution(getOperator().domain());
}


TSFVector TSFMatrixProblem::getKnownSolution(const TSFVectorSpace& space) const 
{
	TSFVector soln = space.createMember();
	fillKnownSolution(soln);
	return soln;
}

void TSFMatrixProblem::fillKnownSolution(TSFVector& soln) const
{
	TSFError::raise("TSFMatrixProblem::fillKnownSolution not implemented");
}

