#include "TSFSumNonlinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"

using namespace TSF;


TSFSumNonlinearOperator::
TSFSumNonlinearOperator(const TSFNonlinearOperator& left, 
												const TSFNonlinearOperator& right, 
												bool subtraction)
	: TSFNonlinearOperatorBase(left.domain(), left.range()),
		left_(left), right_(right), subtraction_(subtraction)
{
	if (left_.domain() != right_.domain())
		{
			TSFError::raise("domain mismatch in TSFSumNonlinearOperator ctor");
		}
	if (left_.range() != right_.range())
		{
			TSFError::raise("range mismatch in TSFSumNonlinearOperator ctor");
		}
}

void TSFSumNonlinearOperator::apply(const TSFVector& in, 
																		TSFVector& out) const
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

TSFLinearOperator TSFSumNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	TSFLinearOperator leftDeriv = left_.derivative(evalPt);
	TSFLinearOperator rightDeriv = right_.derivative(evalPt);

	if (subtraction_)
		{
			return leftDeriv - rightDeriv;
		}
	else
		{
			return leftDeriv + rightDeriv;
		}
	
}

void TSFSumNonlinearOperator::print(ostream& os) const 
{
	os << "(" << left_;
	if (subtraction_)
		{
			os << " - " ;
		}
	else
		{
			os << " + " ;
		}
	os << right_ << ")";
} 


