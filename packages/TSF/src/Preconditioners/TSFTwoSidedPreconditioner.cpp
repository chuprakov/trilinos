#include "TSFTwoSidedPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"

using namespace TSF;
using std::string;

TSFTwoSidedPreconditioner::TSFTwoSidedPreconditioner(const TSFLinearOperator& leftOp,
																						 const TSFLinearOperator& rightOp)
	: TSFPreconditionerBase(), leftOp_(leftOp), rightOp_(rightOp)
{}

string TSFTwoSidedPreconditioner::toString() const
{
	return "TSFTwoSidedPreconditioner(M1^-1 = " + leftOp_.toString() 
		+ ", M2^-1 = " + rightOp_.toString() + ")";
}
