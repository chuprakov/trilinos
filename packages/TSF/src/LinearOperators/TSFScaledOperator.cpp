#include "TSFScaledOperator.h"
#include "TSFSmartPtr.h"
#include "TSFError.h"
#include "TSFUtils.h"

using namespace TSF;


TSFScaledOperator::TSFScaledOperator(const TSFLinearOperator& op, 
																		 const TSFReal& scale)
	: TSFLinearOperatorBase(op.domain(), op.range()), op_(op), scale_(scale)
{}

void TSFScaledOperator::apply(const TSFVector& in, 
															TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
			out.scalarMult(scale_, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFScaledOperator::apply()");
		}
}

void TSFScaledOperator::applyAdjoint(const TSFVector& in, 
																		 TSFVector& out) const
{
	try
		{
			op_.applyAdjoint(in, out);
			out.scalarMult(scale_, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFScaledOperator::applyAdjoint()");
		}
}

void TSFScaledOperator::applyInverse(const TSFVector& in, 
																		 TSFVector& out) const
{
	try
		{
			if (TSFUtils::chop(scale_)==0)
				{
					TSFError::raise("zero divisor in TSFScaledOperator::applyInverse");
				}
			op_.applyInverse(in, out);
			out.scalarMult(1.0/scale_, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFScaledOperator::applyInverse()");
		}
}

