#include "TSFLinearProblem.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFLinearOperatorBase.h"


using namespace TSF;


TSFLinearProblem::TSFLinearProblem(TSFLinearProblemBase* ptr)
	: ptr_(ptr)
{}


TSFVector TSFLinearProblem::getRHS() const 
{
	return ptr_->getRHS();
}

TSFVector TSFLinearProblem::getKnownSolution() const 
{
	return ptr_->getKnownSolution();
}

TSFLinearOperator TSFLinearProblem::getOperator() const 
{
	return ptr_->getOperator();
}


