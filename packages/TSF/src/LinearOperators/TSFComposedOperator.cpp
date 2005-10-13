#include "TSFComposedOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"

using namespace TSF;


TSFComposedOperator::TSFComposedOperator(const TSFLinearOperator& left, 
																				 const TSFLinearOperator& right)
	: TSFLinearOperatorBase(right.domain(), left.range()),
		left_(left), right_(right)
{
	if (left_.domain() != right_.range())
		{
			TSFError::raise("domain-range mismatch in TSFComposedOperator ctor");
		}
}

void TSFComposedOperator::apply(const TSFVector& in, 
																TSFVector& out) const
{
  /* apply operators in order from right to left */
  TSFVector tmp = right_.range().createMember();
  right_.apply(in, tmp);
  left_.apply(tmp, out);
}

void TSFComposedOperator::applyAdjoint(const TSFVector& in, 
																			 TSFVector& out) const
{
	try
		{
			/* adj(A*B) = adj(B)*adj(A). Apply adjoint operators in order from 
			 * left to right */
			TSFVector tmp = left_.domain().createMember();
			left_.applyAdjoint(in, tmp);
			right_.applyAdjoint(tmp, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFComposedOperator::applyAdjoint()");
		}
}

void TSFComposedOperator::applyInverse(const TSFVector& in, 
																			 TSFVector& out) const
{
	try
		{
			/* inv(A*B) = inv(B)*inv(A). Apply inverse operators in order from 
			 * left to right */
			TSFVector tmp = left_.domain().createMember();
			left_.applyInverse(in, tmp);
			right_.applyInverse(tmp, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFComposedOperator::applyInverse()");
		}
}

