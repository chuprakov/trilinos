#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorSpace.h"

#include "TSFVector.h"
#include "DenseSerialVector.h"
#include "TSFError.h"

using namespace TSF;

using std::string;
using std::ostream;

double TSFVectorBase::dummyElement_ = 0.0;

TSFVectorBase::TSFVectorBase(const TSFVectorSpace& space)
	: space_(space)
{;}


TSFVectorBase::~TSFVectorBase()
{;}

int TSFVectorBase::numBlocks() const 
{
	return 1;
}

void TSFVectorBase::getBlock(int i, const TSFVector& self, 
														 TSFVector& sub) const 
{
	if (i==0) 
		{
			sub = self;
		}
	else
		{
			TSFError::raise("TSFVectorBase::getBlock() index error");
		}
}

void TSFVectorBase::setBlock(int /*i*/, const TSFVector& /*sub*/)
{
	TSFError::raise("TSFVectorBase::setBlock()");
}

