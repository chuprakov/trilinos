#include "TSFConstantNonlinearOperator.h"
#include "TSFZeroOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"

using namespace TSF;


TSFConstantNonlinearOperator::
TSFConstantNonlinearOperator(const TSFVectorSpace& domain,
														 const TSFVector& value)
	: TSFNonlinearOperatorBase(domain, value.space()),
		value_(value)
{}

void TSFConstantNonlinearOperator::apply(const TSFVector& /* in */, 
																			 TSFVector& out) const
{
	out = value_.copy();
}

TSFLinearOperator TSFConstantNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	return new TSFZeroOperator(domain(), range());
}

void TSFConstantNonlinearOperator::print(ostream& os) const 
{
	os << "TSFConstantNonlinearOperator(" << value_ << ")";
} 


