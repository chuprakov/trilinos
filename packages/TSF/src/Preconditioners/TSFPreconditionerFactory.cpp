#include "TSFPreconditionerFactory.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"


using namespace TSF;
using std::string;

TSFPreconditionerFactory::TSFPreconditionerFactory()
	: ptr_(0)
{}

TSFPreconditionerFactory
::TSFPreconditionerFactory(TSFPreconditionerFactoryBase* ptr)
	: ptr_(ptr)
{}


TSFPreconditioner TSFPreconditionerFactory::createPreconditioner(const TSFLinearOperator& A) const
{
	if (ptr_==0) return TSFPreconditioner(); 
	return ptr_->createPreconditioner(A);
}

string TSFPreconditionerFactory::toString() const
{
	return ptr_->toString();
}

namespace TSF
{
	ostream& operator<<(ostream& os, const TSFPreconditionerFactory& x)
	{
		return os << x.toString();
	}

	string toString(const TSFPreconditionerFactory& x)
	{
			return x.toString();
		}
}
