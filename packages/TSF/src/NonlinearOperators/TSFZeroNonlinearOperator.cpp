#include "TSFZeroNonlinearOperator.h"
#include "TSFZeroOperator.h"

#include "TSFError.h"

using namespace TSF;


TSFZeroNonlinearOperator::
TSFZeroNonlinearOperator(const TSFVectorSpace& domain,
												 const TSFVectorSpace& range)
	: TSFNonlinearOperatorBase(domain, range)
{}

void TSFZeroNonlinearOperator::apply(const TSFVector& in, 
																		TSFVector& out) const
{
	TSFVector tmp = range().createMember();
	tmp.zero();
}

TSFLinearOperator TSFZeroNonlinearOperator::derivative(const TSFVector& /* evalPt */) const
{
	return new TSFZeroOperator(domain(), range());
}

void TSFZeroNonlinearOperator::print(ostream& os) const 
{
	os << "TSFZeroNonlinearOperator()" ;
} 


