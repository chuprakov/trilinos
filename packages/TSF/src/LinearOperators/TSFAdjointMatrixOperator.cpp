#include "TSFAdjointMatrixOperator.h"
#include "TSFAdjointPreconditioner.h"

namespace TSF {

TSFAdjointMatrixOperator::TSFAdjointMatrixOperator( const TSFSmartPtr<TSFMatrixOperator>& mat_op )
	:mat_op_(mat_op_)
{}

void TSFAdjointMatrixOperator::getILUKPreconditioner(
	int fillLevels, int overlapFill, TSFPreconditioner& rtn) const
{
	TSFPreconditioner prec;
	mat_op_->getILUKPreconditioner(fillLevels,overlapFill,prec);
	rtn = new TSFAdjointPreconditioner(prec);
}

bool TSFAdjointMatrixOperator::isFactored() const
{
	return mat_op_->isFactored();
}

void TSFAdjointMatrixOperator::factor()
{
	mat_op_->factor();
}

void TSFAdjointMatrixOperator::apply(const TSFVector& in, TSFVector& out) const
{
	assert(0);
}

void TSFAdjointMatrixOperator::addToRow(
	int globalRowIndex,
	int nCols,
	const int* globalColumnIndices,
	const TSFReal* a)
{
	assert(0);
}

void TSFAdjointMatrixOperator::zero()
{
	assert(0);
}

} // end namespace TSF
