#include "TSFLeftPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"

using namespace TSF;
using std::string;

TSFLeftPreconditioner::TSFLeftPreconditioner(const TSFLinearOperator& leftOp)
	: TSFPreconditionerBase(), leftOp_(leftOp)
{}

string TSFLeftPreconditioner::toString() const
{
	return "TSFLeftPreconditioner(M1^-1 = " + leftOp_.toString() 
		+ ")";
}
