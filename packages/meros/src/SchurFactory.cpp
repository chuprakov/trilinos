#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFOperatorSourceBase.h"


using namespace SPP;
using namespace TSF;

using std::string;

SchurFactory::SchurFactory()
	: ptr_(0)
{}


SchurFactory
::SchurFactory(SchurFactoryBase* ptr)
	: ptr_(ptr)
{}


TSFLinearOperator SchurFactory::getSchurInvApprox(const TSFOperatorSource& opSrc) const
{
	if (ptr_.get()==0) return TSFLinearOperator(); 
	return ptr_->getSchurInvApprox(opSrc);
}

// TSFLinearOperator SchurFactory::getSchurInvApprox(const KayLoghinRightOperatorSource& opSrc) const
// {
// 	if (ptr_.get()==0) return TSFLinearOperator(); 
// 	return ptr_->getSchurInvApprox(opSrc);
// }


string SchurFactory::toString() const
{
	return ptr_->toString();
}

namespace SPP
{
	ostream& operator<<(ostream& os, const SchurFactory& x)
	{
		return os << x.toString();
	}

	string toString(const SchurFactory& x)
	{
			return x.toString();
		}
}
