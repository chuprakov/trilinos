#include "TSFDeferredCopy.h"
#include "TSFVector.h"
#include "TSFVectorBase.h"
#include "TSFError.h"

#include "TSFOut.h"

using namespace TSF;

TSFDeferredCopy::TSFDeferredCopy(const TSFSmartPtr<TSFVectorBase>& ptr)
	: ptr_(ptr)
{}

void TSFDeferredCopy::createCopy(TSFVector& vec) const
{
	vec = ptr_->deepCopy();
}

void TSFDeferredCopy::copyInto(TSFVector& vec) const
{
	TSFVector tmp;
	tmp.smartPtr() = ptr_;
	vec.acceptCopyOf(tmp);
}
