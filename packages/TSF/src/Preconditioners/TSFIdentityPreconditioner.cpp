#include "TSFIdentityPreconditioner.h"

#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"

using namespace TSF;
using std::string;

TSFIdentityPreconditioner::TSFIdentityPreconditioner()
	: TSFPreconditionerBase()
{}


TSFLinearOperator TSFIdentityPreconditioner::left() const
{
	return new TSFIdentityOperator();
}

string TSFIdentityPreconditioner::toString() const 
{
	return "IdentityPreconditioner";
}
