#include "TSFOperatorSource.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFOperatorSourceBase.h"

using namespace TSF;
using std::string;

TSFOperatorSource::TSFOperatorSource()
	: ptr_(0)
{}

TSFOperatorSource
::TSFOperatorSource(TSFOperatorSourceBase* ptr)
	: ptr_(ptr)
{}


TSFLinearOperator TSFOperatorSource::getOp() const
{
	if (ptr_.get()==0) return TSFLinearOperator(); 
	return ptr_->getOp();
}


string TSFOperatorSource::toString() const
{
	return ptr_->toString();
}

namespace TSF
{
	ostream& operator<<(ostream& os, const TSFOperatorSource& x)
	{
		return os << x.toString();
	}

	string toString(const TSFOperatorSource& x)
	{
			return x.toString();
		}
}
