#include "TSFScaledNonlinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"

using namespace TSF;


TSFScaledNonlinearOperator::
TSFScaledNonlinearOperator(const TSFReal& scale,
													 const TSFNonlinearOperator& op)
	: TSFNonlinearOperatorBase(op.domain(), op.range()),
		op_(op), scale_(scale)
{}

void TSFScaledNonlinearOperator::apply(const TSFVector& in, 
																			 TSFVector& out) const
{
	TSFVector tmp = range().createMember();
	
	op_.apply(in, tmp);
	tmp.scalarMult(scale_, out);
}

TSFLinearOperator TSFScaledNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	return scale_ * op_.derivative(evalPt);
}

void TSFScaledNonlinearOperator::print(ostream& os) const 
{
	os << scale_ << " * " << op_;
} 


