#include "TSFGeneralizedIndex.h"

using namespace TSF;

TSFGeneralizedIndex::TSFGeneralizedIndex(const TSFGeneralizedIndex& index, int i)
	: indices_(index.indices_)
{
	indices_.push(i);
}

TSFGeneralizedIndex::TSFGeneralizedIndex(int i)
	: indices_()
{
	indices_.push(i);
}

TSFGeneralizedIndex::TSFGeneralizedIndex()
	: indices_()
{}

TSFGeneralizedIndex::TSFGeneralizedIndex(const TSFStack<int>& stack)
	: indices_(stack)
{}

TSFGeneralizedIndex TSFGeneralizedIndex::remainder() const 
{
	TSFGeneralizedIndex rtn(indices_);
	rtn.indices_.pop();
	return rtn;
}
