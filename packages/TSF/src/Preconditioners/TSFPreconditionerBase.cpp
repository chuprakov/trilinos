#include "TSFPreconditionerBase.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"

using namespace TSF;
using std::string;

TSFLinearOperator TSFPreconditionerBase::left() const
{
	TSFError::raise("TSFPreconditionerBase::left() undefined for base class");
	return TSFLinearOperator(); // -Wall
}

TSFLinearOperator TSFPreconditionerBase::right() const
{
	TSFError::raise("TSFPreconditionerBase::right() undefined for base class");
	return TSFLinearOperator(); // -Wall
}


