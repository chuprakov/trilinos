#include "TSFPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"

using namespace TSF;
using std::string;

TSFPreconditioner::TSFPreconditioner()
	: ptr_(0)
{}

TSFPreconditioner::TSFPreconditioner(TSFPreconditionerBase* ptr)
	: ptr_(ptr)
{}

bool TSFPreconditioner::hasLeft() const 
{
	if (ptr_.get()==0) return false;
	else return ptr_->hasLeft();
}

bool TSFPreconditioner::hasRight() const 
{
	if (ptr_.get()==0) return false;
	else return ptr_->hasRight();
}




