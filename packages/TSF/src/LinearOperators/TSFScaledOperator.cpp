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


void TSFScaledOperator::getRow(int row, TSFArray<int>& indices, 
                               TSFArray<TSFReal>& values) const
{
  op_.getRow(row, indices, values);
  for (int i = 0; i < indices.size(); i++)
    {
      values[i] = values[i] * scale_;
    }
}

TSFLinearOperator* TSFScaledOperator::getTranspose() 
{
  opTrp_ = new TSFScaledOperator(op_.getTranspose(), scale_);
  return &opTrp_;
}


