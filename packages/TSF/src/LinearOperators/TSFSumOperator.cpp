#include "TSFSumOperator.h"

#include "TSFError.h"

using namespace TSF;


TSFSumOperator::TSFSumOperator(const TSFLinearOperator& left, 
															 const TSFLinearOperator& right, 
															 bool subtraction)
	: TSFLinearOperatorBase(left.domain(), left.range()),
		left_(left), right_(right), subtraction_(subtraction)
{
	if (left_.domain() != right_.domain())
		{
			TSFError::raise("domain mismatch in TSFSumOperator ctor");
		}
	if (left_.range() != right_.range())
		{
			TSFError::raise("range mismatch in TSFSumOperator ctor");
		}
}

void TSFSumOperator::apply(const TSFVector& in, TSFVector& out) const
{
	try
		{
			TSFVector tmp = out.space().createMember();
			
			left_.apply(in, out);
			right_.apply(in, tmp);
			if (subtraction_) 
				{
					out.subtract(out, tmp);
				}
			else
				{
					out.add(out, tmp);
				}
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFSumOperator::apply()");
		}
}
