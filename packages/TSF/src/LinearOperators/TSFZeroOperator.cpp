#include "TSFZeroOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"

using namespace TSF;


TSFZeroOperator::TSFZeroOperator(const TSFVectorSpace& domain,
																 const TSFVectorSpace& range)
	: TSFLinearOperatorBase(domain, range)
{}

void TSFZeroOperator::apply(const TSFVector& /* in */, 
														TSFVector& out) const
{
	try
		{
			out = range().createMember();
			out.zero();
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFIdentityOperator::apply()");
		}
}

void TSFZeroOperator::applyAdjoint(const TSFVector& in, 
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

void TSFZeroOperator::applyInverse(const TSFVector& /* in */, 
																	 TSFVector& /* out */) const
{
	TSFError::raise("Attempted to apply inverse of TSFZeroOperator");
}

void TSFZeroOperator::print(ostream& os) const 
{
	os << "<TSFZeroOperator/>" << endl;
}
