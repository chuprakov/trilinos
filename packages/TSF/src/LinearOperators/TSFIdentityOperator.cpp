#include "TSFIdentityOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"

using namespace TSF;


TSFIdentityOperator::TSFIdentityOperator(const TSFVectorSpace& space)
	: TSFLinearOperatorBase(space, space)
{}

void TSFIdentityOperator::apply(const TSFVector& in, 
																TSFVector& out) const
{
	try
		{
			out = in.copy();
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFIdentityOperator::apply()");
		}
}

void TSFIdentityOperator::applyAdjoint(const TSFVector& in, 
																			 TSFVector& out) const
{
	try
		{
			apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFIdentityOperator::applyAdjoint()");
		}
}

void TSFIdentityOperator::applyInverse(const TSFVector& in, 
																			 TSFVector& out) const
{
	try
		{
			apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFIdentityOperator::applyInverse()");
		}
}

