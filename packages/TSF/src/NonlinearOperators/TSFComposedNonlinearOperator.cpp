#include "TSFComposedNonlinearOperator.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"


using namespace TSF;

TSFComposedNonlinearOperator::TSFComposedNonlinearOperator(const TSFNonlinearOperator& left, 
																													 const TSFNonlinearOperator& right)
	: TSFNonlinearOperatorBase(right.domain(), left.range()),
				left_(left), right_(right) 
{;}

void TSFComposedNonlinearOperator::apply(const TSFVector& arg, TSFVector& out) const
{
	TSFVector tmp;
	right_.apply(arg, tmp);
	left_.apply(tmp, out);
}

TSFLinearOperator TSFComposedNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	TSFVector rightVal;
	right_.apply(evalPt, rightVal);
	return left_.derivative(rightVal) * right_.derivative(evalPt);
}

void TSFComposedNonlinearOperator::print(ostream& os) const 
{
	os << "ComposedNonlinearOperator[left=" << left_ << ", right=" 
		 << right_ << "]";
}








