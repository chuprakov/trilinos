#include "TSFAddVectorNonlinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"

using namespace TSF;


TSFAddVectorNonlinearOperator::
TSFAddVectorNonlinearOperator(const TSFNonlinearOperator& op, 
															const TSFVector& vec, 
															bool subtraction)
	: TSFNonlinearOperatorBase(op.domain(), op.range()),
		op_(op), vec_(vec), subtraction_(subtraction)
{
	if (op_.range() != vec.space())
		{
			TSFError::raise("range mismatch in TSFAddVectorNonlinearOperator ctor");
		}
}

void TSFAddVectorNonlinearOperator::apply(const TSFVector& in, 
																					TSFVector& out) const
{
	TSFVector tmp = out.space().createMember();
	
	op_.apply(in, out);

	if (subtraction_)
		{
			out.subtract(out, vec_);
		}
	else
		{
			out.add(out, vec_);
		}
}

TSFLinearOperator TSFAddVectorNonlinearOperator::derivative(const TSFVector& evalPt) const
{
	return op_.derivative(evalPt);
}

void TSFAddVectorNonlinearOperator::print(ostream& os) const 
{
	os << "TSFAddVectorNonlinearOperator";
} 


