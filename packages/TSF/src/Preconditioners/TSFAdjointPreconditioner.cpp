#include "TSFAdjointPreconditioner.h"
#include "TSFLinearOperatorBase.h"

namespace TSF {

TSFAdjointPreconditioner::TSFAdjointPreconditioner(const TSFPreconditioner& prec)
	:prec_(prec)
{}

TSFLinearOperator TSFAdjointPreconditioner::left() const
{
	return prec_.right().adjoint();
}

TSFLinearOperator TSFAdjointPreconditioner::right() const
{
	return prec_.left().adjoint();
}

bool TSFAdjointPreconditioner::isTwoSided() const
{
	return prec_.isTwoSided();
}

bool TSFAdjointPreconditioner::hasLeft() const
{
	return prec_.hasRight();
}

bool TSFAdjointPreconditioner::hasRight() const
{
	return prec_.hasLeft();
}

bool TSFAdjointPreconditioner::isIdentity() const
{
	return prec_.isIdentity();
}

string TSFAdjointPreconditioner::toString() const
{
	return prec_.toString();
}

} // end namespace TSF
