#include "TSFEmptyVectorSpace.h"
#include "TSFVector.h"
#include "TSFError.h"


using namespace TSF;

TSFVectorSpaceBase* TSFEmptyVectorSpace::deepCopy() const 
{
	return new TSFEmptyVectorSpace();
}

TSFVectorBase* TSFEmptyVectorSpace::createMember(const TSFVectorSpace& /* handle */) const 
{
	TSFError::raise("TSFEmptyVectorSpace::createMember should not be called");
	return 0;
}

