#include "TSFNonlinearOperator.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFLinearNonlinearOperator.h"
#include "TSFSumNonlinearOperator.h"
#include "TSFScaledNonlinearOperator.h"
#include "TSFComposedNonlinearOperator.h"
#include "TSFZeroNonlinearOperator.h"
#include "TSFAddVectorNonlinearOperator.h"

using namespace TSF;

TSFNonlinearOperator::TSFNonlinearOperator()
	: ptr_()
{}


TSFNonlinearOperator::TSFNonlinearOperator(TSFNonlinearOperatorBase* ptr)
	: ptr_(ptr)
{}


TSFNonlinearOperator::TSFNonlinearOperator(const TSFLinearOperator& linOp)
	: ptr_(new TSFLinearNonlinearOperator(linOp))
{}


const TSFVectorSpace& TSFNonlinearOperator::domain() const
{
	return ptr_->domain();
}

const TSFVectorSpace& TSFNonlinearOperator::range() const
{
	return ptr_->range();
}

void TSFNonlinearOperator::apply(const TSFVector& arg, TSFVector& out) const
{
	ptr_->apply(arg, out);
}

TSFLinearOperator TSFNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	return ptr_->derivative(evalPt);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator-() const
{
	return new TSFScaledNonlinearOperator(-1.0, *this);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator*(const TSFReal& scale) const
{
	return new TSFScaledNonlinearOperator(scale, *this);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator+(const TSFNonlinearOperator& op) const
{
	return new TSFSumNonlinearOperator(*this, op, 1);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator-(const TSFNonlinearOperator& op) const
{
	return new TSFSumNonlinearOperator(*this, op, -1);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator+(const TSFVector& v) const
{
	return new TSFAddVectorNonlinearOperator(*this, v, 1);
}

TSFNonlinearOperator 
TSFNonlinearOperator::operator-(const TSFVector& v) const
{
	return new TSFAddVectorNonlinearOperator(*this, v, -1);
}

TSFNonlinearOperator
TSFNonlinearOperator::compose(const TSFNonlinearOperator& other) const 
{
	return new TSFComposedNonlinearOperator(*this, other);
}

void TSFNonlinearOperator::print(ostream& os) const 
{
	ptr_->print(os);
}








