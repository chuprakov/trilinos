#include "TSFRightPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"

using namespace TSF;
using std::string;

TSFRightPreconditioner::TSFRightPreconditioner(const TSFLinearOperator& rightOp)
	: TSFPreconditionerBase(), rightOp_(rightOp)
{}

string TSFRightPreconditioner::toString() const
{
	return "TSFRightPreconditioner(M1^-1 = " + rightOp_.toString() 
		+ ")";
}
